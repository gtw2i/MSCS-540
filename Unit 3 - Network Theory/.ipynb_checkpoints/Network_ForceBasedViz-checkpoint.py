# Network_ForceBasedViz.py
# ------------------------------------------------------------
# Install (once):
#   pip install streamlit pyvis networkx
#
# Run:
#   streamlit run Network_ForceBasedViz.py
# ------------------------------------------------------------

import os
import json
import tempfile
import numpy as np
import streamlit as st
import networkx as nx
from pyvis.network import Network

st.set_page_config(page_title="Interactive Force Graph", layout="wide")

# ============================================================
# Sidebar controls (ALL UI lives here)
# ============================================================
with st.sidebar:
    st.title("Force Graph Controls")

    kind = st.selectbox(
        "Graph type",
        [
            "path",
            "cycle",
            "complete",
            "bipartite (random)",
            "complete bipartite",
            "star",
            "wheel",
            "tree (random)",
            "grid (planar)",
            "planar (random planar)",
            "nonplanar (K3,3)",
            "erdos-renyi",
            "watts-strogatz",
            "barabasi-albert",
        ],
        index=9,
    )

    st.divider()
    st.subheader("Size / Parameters")

    seed = st.number_input("Seed", min_value=0, max_value=10_000_000, value=42, step=1)
    n = st.slider("Number of nodes", 5, 150, 30, 1)

    p = None
    k = None
    m = None
    left_size = None

    if kind == "erdos-renyi":
        p = st.slider("Edge probability p", 0.01, 0.80, 0.10, 0.01)

    if kind in ("bipartite (random)", "complete bipartite"):
        left_size = st.slider("Left partition size", 1, max(2, n - 1), max(2, n // 2), 1)
        if kind == "bipartite (random)":
            p = st.slider("Bipartite edge probability p", 0.01, 0.95, 0.25, 0.01)

    if kind == "watts-strogatz":
        k = st.slider("k (nearest neighbors; even)", 2, min(20, max(2, n - 1)), 6, 1)
        if k % 2 == 1:
            k += 1
        p = st.slider("rewire probability p", 0.0, 1.0, 0.25, 0.01)

    if kind == "barabasi-albert":
        m = st.slider("m (edges per new node)", 1, min(20, max(1, n - 1)), 2, 1)

    st.divider()
    st.subheader("Interaction / Physics")

    show_labels = st.checkbox("Show node labels", value=False)
    physics = st.checkbox("Enable physics", value=True)
    drag_nodes = st.checkbox("Allow dragging nodes", value=True)

    grav = st.slider("Repulsion (gravitationalConstant)", -30000, -500, -8000, 500)
    spring_len = st.slider("Spring length", 20, 300, 120, 5)
    spring_k = st.slider("Spring constant", 1, 100, 4, 1) / 100.0
    damping = st.slider("Damping", 0, 90, 25, 1) / 100.0
    min_vel = st.slider("Min velocity", 0.0, 5.0, 0.75, 0.05)

    st.divider()
    st.subheader("Style")

    node_size = st.slider("Node size", 5, 60, 14, 1)
    edge_width = st.slider("Edge width", 1, 10, 2, 1)
    edge_alpha = st.slider("Edge opacity", 0.05, 1.0, 0.65, 0.05)

    view_height = st.slider("Viewer height (px)", 500, 1100, 820, 10)

    with st.expander("Graph stats (optional)"):
        show_stats = st.checkbox("Show stats here", value=True)

# ============================================================
# Graph construction
# ============================================================
def build_graph(kind: str, n: int, seed: int, p=None, k=None, m=None, left_size=None) -> nx.Graph:
    if kind == "path":
        return nx.path_graph(n)

    if kind == "cycle":
        return nx.cycle_graph(n)

    if kind == "complete":
        return nx.complete_graph(n)

    if kind == "star":
        return nx.star_graph(max(1, n - 1))  # k+1 nodes

    if kind == "wheel":
        return nx.wheel_graph(max(4, n))

    if kind == "tree (random)":
        return nx.random_tree(n, seed=seed)

    if kind == "grid (planar)":
        rows = int(np.floor(np.sqrt(n)))
        cols = int(np.ceil(n / max(1, rows)))
        G = nx.grid_2d_graph(rows, cols)
        G = nx.convert_node_labels_to_integers(G, ordering="sorted")
        if G.number_of_nodes() > n:
            keep = list(range(n))
            return G.subgraph(keep).copy()
        return G

    if kind == "planar (random planar)":
        try:
            return nx.random_planar_graph(n, seed=seed)
        except Exception:
            rows = int(np.floor(np.sqrt(n)))
            cols = int(np.ceil(n / max(1, rows)))
            Gt = nx.triangular_lattice_graph(rows, cols)
            Gt = nx.convert_node_labels_to_integers(Gt, ordering="sorted")
            if Gt.number_of_nodes() > n:
                keep = list(range(n))
                return Gt.subgraph(keep).copy()
            return Gt

    if kind == "nonplanar (K3,3)":
        return nx.complete_bipartite_graph(3, 3)

    if kind == "erdos-renyi":
        return nx.erdos_renyi_graph(n=n, p=float(p), seed=seed)

    if kind == "watts-strogatz":
        kk = int(k)
        kk = min(kk, max(2, n - 1))
        if kk % 2 == 1:
            kk -= 1
        return nx.watts_strogatz_graph(n=n, k=kk, p=float(p), seed=seed)

    if kind == "barabasi-albert":
        mm = int(m)
        mm = min(mm, max(1, n - 1))
        return nx.barabasi_albert_graph(n=n, m=mm, seed=seed)

    if kind == "bipartite (random)":
        a = int(left_size)
        b = max(1, n - a)
        return nx.bipartite.random_graph(a, b, float(p), seed=seed)

    if kind == "complete bipartite":
        a = int(left_size)
        b = max(1, n - a)
        return nx.complete_bipartite_graph(a, b)

    raise ValueError(f"Unknown kind: {kind}")

def initial_pos(kind: str, G: nx.Graph, seed: int, left_size=None) -> dict:
    if kind == "path":
        nodes = sorted(G.nodes())
        return {u: (i * 60.0, 0.0) for i, u in enumerate(nodes)}

    if kind in ("cycle", "complete"):
        pos = nx.circular_layout(G, scale=220.0)
        return {u: (float(pos[u][0]), float(pos[u][1])) for u in G.nodes()}

    if kind == "star":
        center = 0
        leaves = [u for u in G.nodes() if u != center]
        pos = {center: (0.0, 0.0)}
        angles = np.linspace(0, 2*np.pi, len(leaves), endpoint=False)
        for a, u in zip(angles, leaves):
            pos[u] = (220.0 * np.cos(a), 220.0 * np.sin(a))
        return pos

    if kind == "wheel":
        hub = 0
        rim = [u for u in G.nodes() if u != hub]
        pos = {hub: (0.0, 0.0)}
        angles = np.linspace(0, 2*np.pi, len(rim), endpoint=False)
        for a, u in zip(angles, rim):
            pos[u] = (260.0 * np.cos(a), 260.0 * np.sin(a))
        return pos

    if kind == "grid (planar)":
        N = G.number_of_nodes()
        rows = int(np.floor(np.sqrt(N)))
        cols = int(np.ceil(N / max(1, rows)))
        pos = {}
        for u in G.nodes():
            r = u // cols
            c = u % cols
            pos[u] = (c * 80.0, -r * 80.0)
        return pos

    if kind in ("bipartite (random)", "complete bipartite", "nonplanar (K3,3)"):
        a = int(left_size) if left_size is not None else 3
        left = list(range(a))
        pos = nx.bipartite_layout(G, left, scale=260.0)
        return {u: (float(pos[u][0]), float(pos[u][1])) for u in G.nodes()}

    if kind == "tree (random)":
        root = 0
        DG = nx.bfs_tree(G, source=root)
        layers = {}
        for u in nx.topological_sort(DG):
            d = nx.shortest_path_length(DG, root, u)
            layers.setdefault(d, []).append(u)

        pos = {}
        for d, nodes in layers.items():
            k = len(nodes)
            xs = (np.arange(k) - (k - 1) / 2.0) * 80.0
            y = -d * 90.0
            for x, u in zip(xs, nodes):
                pos[u] = (float(x), float(y))
        return pos

    if kind == "planar (random planar)":
        try:
            pos = nx.planar_layout(G, scale=260.0)
            return {u: (float(pos[u][0]), float(pos[u][1])) for u in G.nodes()}
        except Exception:
            pass

    pos = nx.spring_layout(G, seed=seed, scale=260.0)
    return {u: (float(pos[u][0]), float(pos[u][1])) for u in G.nodes()}

# Build graph
G = build_graph(kind, n, int(seed), p=p, k=k, m=m, left_size=left_size)

# Connect components lightly for nicer physics (if needed)
UG = G.to_undirected()
if UG.number_of_nodes() > 0 and not nx.is_connected(UG):
    comps = list(nx.connected_components(UG))
    for i in range(len(comps) - 1):
        a = next(iter(comps[i]))
        b = next(iter(comps[i + 1]))
        G.add_edge(a, b)
    UG = G.to_undirected()

pos0 = initial_pos(kind, G, int(seed), left_size=left_size)

# Optional stats in sidebar only
with st.sidebar:
    if "show_stats" in locals() and show_stats:
        st.write(
            {
                "nodes": int(G.number_of_nodes()),
                "edges": int(G.number_of_edges()),
                "connected (undirected)": bool(nx.is_connected(UG)) if UG.number_of_nodes() > 0 else True,
                "density": float(nx.density(UG)) if UG.number_of_nodes() > 1 else 0.0,
            }
        )
        st.caption("Drag nodes in the main view. With physics on, releasing a node re-stabilizes the layout.")

# ============================================================
# PyVis network
# ============================================================
net = Network(height=f"{int(view_height)}px", width="100%", bgcolor="#111111", font_color="white", directed=False)
net.barnes_hut()

for u in G.nodes():
    x, y = pos0.get(u, (0.0, 0.0))
    net.add_node(
        u,
        label=str(u) if show_labels else "",
        title=f"node {u}",
        x=x,
        y=y,
        physics=True,
    )

for u, v in G.edges():
    net.add_edge(u, v)

options = {
    "interaction": {
        "dragNodes": bool(drag_nodes),
        "hover": True,
        "multiselect": True,
        "navigationButtons": True,
        "keyboard": True,
    },
    "physics": {
        "enabled": bool(physics),
        "barnesHut": {
            "gravitationalConstant": int(grav),
            "springLength": int(spring_len),
            "springConstant": float(spring_k),
            "damping": float(damping),
            "avoidOverlap": 0.2,
        },
        "minVelocity": float(min_vel),
    },
    "nodes": {
        "shape": "dot",
        "size": int(node_size),
        "borderWidth": 2,
        "color": {
            "background": "red",
            "border": "white",
            "highlight": {"background": "red", "border": "#ffd166"},
        },
    },
    "edges": {
        "width": int(edge_width),
        "color": {"color": f"rgba(255,255,255,{float(edge_alpha)})"},
        "smooth": False,
    },
}
net.set_options(json.dumps(options))

# ============================================================
# Main area: ONLY the graph visualization
# ============================================================
with tempfile.TemporaryDirectory() as tmpdir:
    html_path = os.path.join(tmpdir, "graph.html")
    net.write_html(html_path, open_browser=False)
    html = open(html_path, "r", encoding="utf-8").read()

st.components.v1.html(html, height=int(view_height) + 30, scrolling=True)