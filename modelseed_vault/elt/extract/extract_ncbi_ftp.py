

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
