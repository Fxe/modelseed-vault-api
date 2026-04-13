import requests
import json
from modelseed_vault.core.transform_graph import Node, HashNode


class Vault:

    def __init__(self, url="http://192.168.1.22:12022/"):
        self.session = requests.Session()
        self.url = url

    def close(self):
        self.session.close()

    def cobra_get_model(self, model_id: str):
        url = f"{self.url}/cobra/model/{model_id}"
        res = self.session.get(url)
        if res.status_code != 200:
            return None
        if res.content:
            return json.loads(res.content)
        return None

    def register(self, node_type):
        url = f"{self.url}/graph/node/constraint"
        params = {
            "type": node_type
        }
        headers = {
            "accept": "*/*"
        }
        data = ""  # empty body

        return requests.post(url, headers=headers, params=params, data=data)

    def get_constraints(self):
        url = f"{self.url}/graph/node/constraint"
        res = requests.get(url)
        if res.status_code != 200:
            raise
        return json.loads(res.content)

    def add_node(self, node_type, node_id, properties):
        url = f"{self.url}/graph/node/{node_type}/{node_id}"
        headers = {
            "accept": "*/*",
            "Content-Type": "application/json",
        }
        return self.session.post(url, headers=headers, json=properties)

    """
    def add_node(self, node_type, node_id, properties):
        url = f"{self.url}/graph/node/{node_type}/{node_id}"
        res = self.session.post(
            url,
            headers=self.headers,
            params=properties,
            data="",
            timeout=(5, 30),
        )
        res.raise_for_status()
        return res
    """

    def add_node2(self, node: Node):
        node_type = node.primary_label
        node_id = node.key
        properties = node.data
        """
        example
        curl -X 'POST' \
  'http://localhost:18080/graph/node/Fruit/banana?labels=Yellow&labels=Tree' \
  -H 'accept: */*' \
  -H 'Content-Type: application/json' \
  -d '{
  "additionalProp1": "x",
  "additionalProp2": "y",
  "additionalProp3": "z"
}'
        add node.labels
        """
        # node.labels

        url = f"{self.url}/graph/node/{node_type}/{node_id}"
        headers = {
            "accept": "*/*",
            "Content-Type": "application/json",
        }
        return self.session.post(url, headers=headers, json=properties)

    def bulk_add_nodes(self, nodes):
        url = f"{self.url}/graph/bulk/nodes"
        headers = {
            "accept": "*/*",
            "Content-Type": "application/json",
        }
        data = []
        for node in nodes:
            node_labels = []
            if node.labels is not None:
                node_labels = [x for x in node.labels]
            data.append({
                "type": node.primary_label,
                "id": node.key,
                "labels": node_labels,
                "properties": node.data
            })

        response = self.session.post(url, headers=headers, json=data)
        response.raise_for_status()

        return response.json()

    def get_node(self, node_type, node_id) -> Node | None:
        url = f"{self.url}/graph/node/{node_type}/{node_id}"
        res = self.session.get(url)
        res.raise_for_status()

        if res.content:
            return json.loads(res.content)
        return None

    def list_nodes(self, node_type: str) -> Node | None:
        url = f"{self.url}/graph/node/{node_type}"
        res = requests.get(url)
        res.raise_for_status()

        return json.loads(res.content)

    def add_protein(self, protein_node: HashNode) -> str:
        url = f"{self.url}/protein/"
        headers = {
            "accept": "*/*",
            "Content-Type": "application/json",
        }
        res = self.session.post(url, headers=headers, json=protein_node._value)
        res.raise_for_status()
        return res.content.decode('utf-8')

    def get_protein_by_sha256(self, sha256: str):
        url = f"{self.url}/protein/sha256/{sha256}"
        res = self.session.get(url)
        if res.status_code != 200:
            return None
        if res.content:
            return json.loads(res.content)
        return None

    def add_edge(self, node_from_eid, node_to_eid, edge_type, props=None):
        from urllib.parse import quote
        src = quote(str(node_from_eid), safe="")
        dst = quote(str(node_to_eid), safe="")
        etype = quote(str(edge_type), safe="")
        url = f"{self.url}/graph/edge/{src}/{dst}/{etype}"
        headers = {
            "accept": "*/*",
            "Content-Type": "application/json",
        }
        return self.session.post(url, headers=headers, json=props)

    def get_node_child(self, node: Node, rel_type=None):
        url = f"{self.url}/graph/node/{node.primary_label}/{node.key}/child"
        if rel_type:
            url += f"?edgeType={rel_type}"
        res = self.session.get(url)
        if res.status_code != 200:
            raise
        return json.loads(res.content)

    def get_node_parent(self, node: Node, rel_type=None):
        url = f"{self.url}/graph/node/{node.primary_label}/{node.key}/parent"
        if rel_type:
            url += f"?edgeType={rel_type}"
        res = self.session.get(url)
        if res.status_code != 200:
            raise
        return json.loads(res.content)
