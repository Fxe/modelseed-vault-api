import re
from ftplib import error_perm, error_proto, error_temp
from io import BytesIO
from minio.error import S3Error


def fetch_ncbi_ftp_data(self, ftp):
    """

    :param ftp: FTP client
    :return:
    """
    ftp_path = self.cwd_ftp_path_rs
    if ftp_path is None:
        ftp_path = self.cwd_ftp_path_gb

    if ftp_path:
        write_path = f'{self.cache_folder}/{ftp_path}'
        os.makedirs(write_path, exist_ok=True)

        ftp.cwd(ftp_path)
        files = ftp.nlst()
        for f in files:
            target_file = f'{write_path}/{f}'
            if f.endswith('_assembly_structure'):  # TODO: implement fetch _assembly_structure
                os.makedirs(target_file, exist_ok=True)
            else:
                with open(f'{write_path}/{f}', 'wb') as fh:
                    ftp.retrbinary(f"RETR {f}", fh.write)


class ExtractNCBIFTP:

    ACCESSION_PATTERN = re.compile(r'^(GCA|GCF)_(\d{9})(\.\d+)?$')
    NCBI_GENOMES_ALL = '/genomes/all'

    def __init__(self):
        pass

    def build_ftp_path(self, ncbi_genome_accession_id: str) -> str:
        """
        Build the NCBI FTP directory holding the assembly folders for an accession.

        GCF_000328665.1 -> /genomes/all/GCF/000/328/665

        :param ncbi_genome_accession_id: NCBI accession_id like GCF_000328665.1
        :return: absolute FTP path of the parent directory
        """
        match = self.ACCESSION_PATTERN.match(ncbi_genome_accession_id.strip())
        if match is None:
            raise ValueError(f'not a valid NCBI assembly accession: {ncbi_genome_accession_id!r}')

        prefix, digits, _ = match.groups()
        return f'{self.NCBI_GENOMES_ALL}/{prefix}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}'

    def resolve_assembly_dir(self, ftp, ncbi_genome_accession_id: str) -> str:
        """
        Resolve the versioned assembly directory, whose name embeds the assembly
        name we do not know up front (GCF_000328665.1_ASM32866v1).

        :param ftp: FTP client
        :param ncbi_genome_accession_id: NCBI accession_id like GCF_000328665.1
        :return: absolute FTP path of the assembly directory
        """
        parent = self.build_ftp_path(ncbi_genome_accession_id)
        accession = ncbi_genome_accession_id.strip()

        try:
            entries = ftp.nlst(parent)
        except (error_perm, error_temp) as err:
            # NCBI answers 450 for a missing directory, other servers use 550
            raise FileNotFoundError(f'no NCBI FTP directory for {accession} at {parent}') from err

        # nlst may return bare names or full paths depending on the server
        names = [e.rsplit('/', 1)[-1] for e in entries]
        matches = sorted(n for n in names if n.startswith(f'{accession}_'))

        if not matches:
            raise FileNotFoundError(f'no assembly directory for {accession} under {parent}')

        # an accession pins one assembly; if the server ever lists more, take the last
        return f'{parent}/{matches[-1]}'

    def list_ftp_dir(self, ftp, path):
        """
        List a remote directory, separating files from sub directories.

        :param ftp: FTP client
        :param path: absolute FTP path
        :return: (files, dirs) as lists of (name, size) and names
        """
        files, dirs = [], []

        try:
            symlinks = []
            for name, facts in ftp.mlsd(path):
                if name in ('.', '..'):
                    continue
                entry_type = facts.get('type')
                if entry_type == 'dir':
                    dirs.append(name)
                elif entry_type == 'file':
                    files.append((name, int(facts['size']) if 'size' in facts else None))
                elif entry_type and entry_type.endswith('=symlink'):
                    # the advertised size is the link target string, not the payload
                    symlinks.append(name)

            # SIZE follows the link: it succeeds for a file, fails for a directory.
            # Symlinked directories are skipped, they duplicate trees carried
            # elsewhere and can point back up and loop.
            if symlinks:
                ftp.voidcmd('TYPE I')  # mlsd left the connection in ASCII
                for name in symlinks:
                    try:
                        files.append((name, ftp.size(f'{path}/{name}')))
                    except error_perm:
                        pass

            return files, dirs
        except (error_perm, error_proto):
            # server without MLSD support, fall back to nlst + SIZE
            files, dirs = [], []  # drop anything MLSD emitted before failing

        entries = ftp.nlst(path)
        ftp.voidcmd('TYPE I')  # nlst left the connection in ASCII, where SIZE is refused

        # LIST exposes the mode bits, the only way to spot symlinks without MLSD
        linked = set()
        try:
            lines = []
            ftp.retrlines(f'LIST {path}', lines.append)
            for line in lines:
                if line.startswith('l'):
                    entry = line.split(None, 8)[-1]
                    linked.add(entry.split(' -> ', 1)[0].rsplit('/', 1)[-1])
            ftp.voidcmd('TYPE I')
        except (error_perm, error_proto, IndexError):
            pass  # no LIST, fall through and treat every entry as a real one

        for entry in entries:
            name = entry.rsplit('/', 1)[-1]
            if name in ('.', '..'):
                continue
            try:
                files.append((name, ftp.size(f'{path}/{name}')))
            except error_perm:
                if name not in linked:  # skip symlinked dirs, see the MLSD branch
                    dirs.append(name)

        return files, dirs

    def stream_ftp_to_minio(self, ftp, minio_client, bucket, remote_file, size, object_name):
        """
        Pipe one FTP file straight into MinIO without staging it on local disk.

        :param ftp: FTP client
        :param minio_client: MinIO client
        :param bucket: target bucket
        :param remote_file: absolute FTP path of the file
        :param size: file size in bytes, or None when the server would not report it
        :param object_name: target object name
        """
        # listing commands (NLST/MLSD) leave the connection in ASCII mode, which would
        # corrupt binary payloads, so force binary again for every transfer
        ftp.voidcmd('TYPE I')

        if size is None:
            # length is required for a single-shot upload, so buffer this one
            buffer = BytesIO()
            ftp.retrbinary(f'RETR {remote_file}', buffer.write)
            size = buffer.tell()
            buffer.seek(0)
            minio_client.put_object(bucket, object_name, buffer, size)
            return

        conn = ftp.transfercmd(f'RETR {remote_file}')
        try:
            minio_client.put_object(bucket, object_name, conn.makefile('rb'), size)
        finally:
            conn.close()
        ftp.voidresp()

    def transfer_genome(self, ncbi_genome_accession_id: str, ftp_client, minio_client, minio_path,
                        skip_existing: bool = True):
        """
        Transfer NCBI to MinIO

        :param ncbi_genome_accession_id: NCBI accession_id like GCF_000328665.1
        :param ftp_client: FTP client already connected to ftp.ncbi.nlm.nih.gov
        :param minio_client: MinIO client
        :param minio_path: target as 'bucket' or 'bucket/prefix'
        :param skip_existing: leave objects already present at the same size alone
        :return: list of object names copied
        """
        bucket, _, prefix = minio_path.strip('/').partition('/')
        if not minio_client.bucket_exists(bucket):
            raise ValueError(f'MinIO bucket does not exist: {bucket}')

        assembly_dir = self.resolve_assembly_dir(ftp_client, ncbi_genome_accession_id)
        assembly_name = assembly_dir.rsplit('/', 1)[-1]

        base = f'{prefix}/{assembly_name}' if prefix else assembly_name
        copied = []

        # walk the assembly directory, _assembly_structure and friends included
        pending = [(assembly_dir, base)]
        while pending:
            remote_dir, object_dir = pending.pop()
            files, dirs = self.list_ftp_dir(ftp_client, remote_dir)

            for name in dirs:
                pending.append((f'{remote_dir}/{name}', f'{object_dir}/{name}'))

            for name, size in files:
                object_name = f'{object_dir}/{name}'

                if skip_existing and size is not None:
                    try:
                        if minio_client.stat_object(bucket, object_name).size == size:
                            continue
                    except S3Error as err:
                        if err.code != 'NoSuchKey':
                            raise

                self.stream_ftp_to_minio(ftp_client, minio_client, bucket,
                                    f'{remote_dir}/{name}', size, object_name)
                copied.append(object_name)

        return copied
