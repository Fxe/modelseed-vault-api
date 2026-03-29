import requests
import json


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

    def add_node2(self, node):
        node_type = node.label
        node_id = node.key
        properties = node.data

        url = f"{self.url}/graph/node/{node_type}/{node_id}"
        headers = {
            "accept": "*/*",
            "Content-Type": "application/json",
        }
        return self.session.post(url, headers=headers, json=properties)

    def get_node(self, node_type, node_id):
        url = f"{self.url}/graph/node/{node_type}/{node_id}"
        res = self.session.get(url)
        if res.status_code != 200:
            return None
        if res.content:
            return json.loads(res.content)
        return None

    def list_nodes(self, node_type):
        url = f"{self.url}/graph/node/{node_type}"
        res = requests.get(url)
        if res.status_code != 200:
            return None
        return json.loads(res.content)

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

    def get_node_child(self, node, rel_type=None):
        url = f"{self.url}/graph/node/{node.label}/{node.key}/child"
        if rel_type:
            url += f"?edgeType={rel_type}"
        res = requests.get(url)
        if res.status_code != 200:
            raise
        return json.loads(res.content)
