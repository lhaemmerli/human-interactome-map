"""Build the human protein interaction map dataset.

Inputs (data/raw/):
  string_physical.txt.gz  STRING v12 physical interactions (human, 9606)
  string_info.txt.gz      STRING protein names / annotations
  uniprot_human.tsv.gz    UniProt reviewed human proteins (keywords + STRING xref)
  proteinatlas.tsv        Human Protein Atlas (protein class, tissue specificity)

Output: docs/data.json
"""

import csv
import gzip
import json
import math
import random
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import networkx as nx

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "docs" / "data.json"

MIN_SCORE = int(sys.argv[1]) if len(sys.argv) > 1 else 700  # STRING high confidence
csv.field_size_limit(10**9)

# ---------------------------------------------------------------- functions
# Each function class is matched from UniProt keywords and/or HPA protein classes.
FUNCTIONS = {
    "Kinase": {"kw": {"Kinase"}},
    "Phosphatase": {"kw": {"Protein phosphatase"}},
    "Transcription regulation": {"kw": {"Transcription regulation", "Activator", "Repressor"},
                                 "hpa": {"Transcription factors"}},
    "Chromatin & epigenetics": {"kw": {"Chromatin regulator"}},
    "DNA replication & repair": {"kw": {"DNA repair", "DNA replication", "DNA damage"}},
    "RNA processing & splicing": {"kw": {"mRNA splicing", "mRNA processing", "rRNA processing",
                                         "tRNA processing", "mRNA transport"}},
    "Translation & ribosome": {"kw": {"Ribosomal protein", "Protein biosynthesis", "Initiation factor",
                                      "Elongation factor"}, "hpa": {"Ribosomal proteins"}},
    "Protein degradation (ubiquitin/protease)": {"kw": {"Ubl conjugation pathway", "Protease"}},
    "Cell cycle & division": {"kw": {"Cell cycle", "Cell division", "Mitosis"}},
    "Apoptosis": {"kw": {"Apoptosis"}},
    "Receptor": {"kw": {"Receptor"}, "hpa": {"G-protein coupled receptors", "Nuclear receptors"}},
    "Ion channel": {"kw": {"Ion channel"}, "hpa": {"Voltage-gated ion channels"}},
    "Transporter": {"hpa": {"Transporters"}},
    "Metabolic enzyme": {"hpa": {"Metabolic proteins"}},
    "Signal transduction (GTPases)": {"kw": {"GTP-binding"}},
    "Hormones, cytokines & growth factors": {"kw": {"Hormone", "Cytokine", "Growth factor"}},
    "Immune response": {"kw": {"Immunity", "Innate immunity", "Adaptive immunity",
                               "Inflammatory response"}},
    "Cytoskeleton & motility": {"kw": {"Cytoskeleton", "Actin-binding", "Microtubule", "Motor protein"}},
    "Cell adhesion & ECM": {"kw": {"Cell adhesion", "Extracellular matrix", "Collagen"}},
    "Vesicle trafficking": {"kw": {"Protein transport", "Endocytosis", "Exocytosis", "ER-Golgi transport"}},
    "Development & differentiation": {"kw": {"Developmental protein", "Differentiation"}},
    "Neurotransmission & synapse": {"kw": {"Neurotransmitter degradation", "Synapse", "Neurogenesis"}},
}
for spec in FUNCTIONS.values():
    spec.setdefault("kw", set())
    spec.setdefault("hpa", set())

# ----------------------------------------------------------------- diseases
# Grouped from the HPA "Disease involvement" column (UniProt disease keywords).
DISEASES = {
    "Cancer": {"Cancer-related genes", "Proto-oncogene", "Tumor suppressor", "Oncogene",
               "Hereditary nonpolyposis colorectal cancer", "Li-Fraumeni syndrome", "Xeroderma pigmentosum"},
    "Alzheimer disease": {"Alzheimer disease"},
    "Parkinson disease": {"Parkinson disease", "Parkinsonism"},
    "Amyotrophic lateral sclerosis (ALS)": {"Amyotrophic lateral sclerosis"},
    "Other neurodegeneration": {"Neurodegeneration", "Spinocerebellar ataxia", "Hereditary spastic paraplegia",
                                "Leukodystrophy", "Neuronal ceroid lipofuscinosis", "Metachromatic leukodystrophy",
                                "Pontocerebellar hypoplasia"},
    "Epilepsy": {"Epilepsy"},
    "Intellectual disability": {"Intellectual disability", "Primary microcephaly", "Lissencephaly"},
    "Autism": {"Autism", "Autism spectrum disorder", "Asperger syndrome"},
    "Schizophrenia": {"Schizophrenia"},
    "Diabetes": {"Diabetes mellitus"},
    "Obesity & lipid disorders": {"Obesity", "Hyperlipidemia", "Congenital generalized lipodystrophy"},
    "Heart & vascular disease": {"Cardiomyopathy", "Atrial fibrillation", "Long QT syndrome", "Short QT syndrome",
                                 "Brugada syndrome", "Aortic aneurysm", "Atherosclerosis", "Atrial septal defect",
                                 "Thrombophilia", "Heterotaxy"},
    "Blood disorders (anemia, bleeding)": {"Hereditary hemolytic anemia", "Fanconi anemia", "Diamond-Blackfan anemia",
                                           "Congenital dyserythropoietic anemia", "Elliptocytosis", "Pyropoikilocytosis",
                                           "Hemophilia", "von Willebrand disease", "Bernard Soulier syndrome",
                                           "Congenital erythrocytosis", "Hemolytic uremic syndrome",
                                           "Hermansky-Pudlak syndrome", "Dyskeratosis congenita"},
    "Immune & autoimmune disease": {"SCID", "Systemic lupus erythematosus", "Chronic granulomatous disease",
                                    "Familial hemophagocytic lymphohistiocytosis", "Aicardi-Goutieres syndrome",
                                    "Asthma", "Allergen", "Autoimmune encephalomyelitis"},
    "Muscle & nerve disease": {"Neuropathy", "Charcot-Marie-Tooth disease", "Dejerine-Sottas syndrome",
                               "Limb-girdle muscular dystrophy", "Congenital muscular dystrophy",
                               "Emery-Dreifuss muscular dystrophy", "Dystroglycanopathy", "Myofibrillar myopathy",
                               "Nemaline myopathy", "Desmin-related myopathy", "Congenital myasthenic syndrome",
                               "Dystonia"},
    "Eye disease & blindness": {"Retinitis pigmentosa", "Cataract", "Cone-rod dystrophy", "Leber congenital amaurosis",
                                "Age-related macular degeneration", "Corneal dystrophy", "Stargardt disease",
                                "Microphthalmia", "Glaucoma", "Congenital stationary night blindness",
                                "Leber hereditary optic neuropathy", "Usher syndrome", "Peters anomaly"},
    "Deafness": {"Deafness", "Non-syndromic deafness", "Usher syndrome", "Waardenburg syndrome"},
    "Mitochondrial disease": {"Primary mitochondrial disease", "Leigh syndrome", "MELAS syndrome",
                              "Progressive external ophthalmoplegia", "Leber hereditary optic neuropathy"},
    "Inborn errors of metabolism": {"Glycogen storage disease", "Congenital disorder of glycosylation",
                                    "Peroxisome biogenesis disorder", "Zellweger syndrome", "Mucopolysaccharidosis",
                                    "Mucolipidosis", "Gangliosidosis", "Gaucher disease", "Niemann-Pick disease",
                                    "Phenylketonuria", "Maple syrup urine disease", "Glutaricaciduria", "Cystinuria",
                                    "Gout", "Rhizomelic chondrodysplasia punctata"},
    "Kidney disease": {"Nephronophthisis", "Alport syndrome", "Bartter syndrome", "Primary hypomagnesemia",
                       "Senior-Loken syndrome", "Hemolytic uremic syndrome", "Diabetes insipidus"},
    "Ciliopathy": {"Ciliopathy", "Primary ciliary dyskinesia", "Kartagener syndrome", "Joubert syndrome",
                   "Bardet-Biedl syndrome", "Meckel syndrome"},
    "Bone & growth disorders": {"Dwarfism", "Osteogenesis imperfecta", "Osteoporosis", "Osteopetrosis",
                                "Craniosynostosis", "Ehlers-Danlos syndrome", "Stickler syndrome",
                                "Hereditary multiple exostoses", "Amelogenesis imperfecta"},
    "Skin & hair disease": {"Ichthyosis", "Palmoplantar keratoderma", "Ectodermal dysplasia", "Epidermolysis bullosa",
                            "Hypotrichosis", "Albinism"},
    "Hormonal disorders": {"Congenital hypothyroidism", "Hypogonadotropic hypogonadism", "Kallmann syndrome",
                           "Cushing syndrome", "Congenital adrenal hyperplasia", "Premature ovarian failure",
                           "Pseudohermaphroditism"},
    "Any inherited disease": {"Disease variant"},
    "Target of FDA-approved drugs": {"FDA approved drug targets"},
}
HIDDEN_DISEASE_TERMS = {"Disease variant"}

# ------------------------------------------------------------- body systems
TISSUE_TO_SYSTEM = {
    "brain": "Nervous", "retina": "Nervous", "choroid plexus": "Nervous",
    "pituitary gland": "Endocrine", "thyroid gland": "Endocrine", "parathyroid gland": "Endocrine",
    "adrenal gland": "Endocrine",
    "esophagus": "Digestive", "stomach 1": "Digestive", "intestine": "Digestive", "liver": "Digestive",
    "gallbladder": "Digestive", "pancreas": "Digestive", "salivary gland": "Digestive", "tongue": "Digestive",
    "heart muscle": "Cardiovascular", "blood vessel": "Cardiovascular",
    "lung": "Respiratory",
    "kidney": "Urinary", "urinary bladder": "Urinary",
    "bone marrow": "Immune & blood", "lymphoid tissue": "Immune & blood",
    "skeletal muscle": "Musculoskeletal", "smooth muscle": "Musculoskeletal",
    "skin 1": "Skin & adipose", "adipose tissue": "Skin & adipose",
    "breast": "Reproductive", "cervix": "Reproductive", "endometrium 1": "Reproductive",
    "epididymis": "Reproductive", "fallopian tube": "Reproductive", "ovary": "Reproductive",
    "placenta": "Reproductive", "prostate": "Reproductive", "seminal vesicle": "Reproductive",
    "testis": "Reproductive", "vagina": "Reproductive",
}
UBIQUITOUS = "Ubiquitous (all tissues)"
SPECIFIC_LEVELS = {"Tissue enriched", "Group enriched", "Tissue enhanced"}


def parse_protein_names(text):
    """'Actin, cytoplasmic 1 (EC 3.6.4.-) (Beta-actin) [Cleaved into: ...]' -> ['Actin, cytoplasmic 1', 'Beta-actin']"""
    text = re.split(r"\s\[(?:Cleaved into|Includes):", text)[0]  # drop sub-chain sections
    out, depth, bdepth, cur = [], 0, 0, ""
    for ch in text:
        if ch in "[]":  # brackets are part of a name: "[Pyruvate dehydrogenase (acetyl-transferring)] kinase"
            bdepth += 1 if ch == "[" else -1
        elif bdepth > 0:
            pass
        elif ch == "(":
            if depth == 0 and cur.strip():
                out.append(cur.strip())
                cur = ""
            depth += 1
            if depth == 1:
                continue
        elif ch == ")":
            depth -= 1
            if depth == 0:
                out.append(cur.strip())
                cur = ""
                continue
        cur += ch
    if cur.strip() and depth == 0:
        out.append(cur.strip())
    return [n for n in dict.fromkeys(out) if n and not n.startswith("EC ")]


def read_tsv(path, opener=open):
    with opener(path, "rt", encoding="utf-8") as fh:
        yield from csv.DictReader(fh, delimiter="\t")


def main():
    # --- STRING names
    info = {}
    for row in read_tsv(RAW / "string_info.txt.gz", gzip.open):
        info[row["#string_protein_id"]] = (row["preferred_name"], row["annotation"])

    # --- UniProt: STRING id -> record
    uniprot, uniprot_by_gene = {}, {}
    for row in read_tsv(RAW / "uniprot_human.tsv.gz", gzip.open):
        for sid in filter(None, row["STRING"].split(";")):
            uniprot.setdefault(sid, row)
        if row["Gene Names (primary)"]:
            uniprot_by_gene.setdefault(row["Gene Names (primary)"], row)

    # --- Pfam domain families per UniProt accession
    pfam_info = {}
    with gzip.open(RAW / "pfam_clans.tsv.gz", "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            pfam_info[f[0]] = (f[3], f[4])  # short name, description
    pfam_of = {}
    for row in read_tsv(RAW / "uniprot_pfam.tsv.gz", gzip.open):
        pfam_of[row["Entry"]] = [x for x in row["Pfam"].split(";") if x in pfam_info]

    # --- HPA: keyed by UniProt accession and gene name
    hpa_by_acc, hpa_by_gene = {}, {}
    for row in read_tsv(RAW / "proteinatlas.tsv"):
        for acc in filter(None, (a.strip() for a in row["Uniprot"].split(","))):
            hpa_by_acc.setdefault(acc, row)
        hpa_by_gene.setdefault(row["Gene"], row)

    # --- Interaction graph
    G = nx.Graph()
    with gzip.open(RAW / "string_physical.txt.gz", "rt") as fh:
        next(fh)  # header: protein1 protein2 combined_score (space-separated)
        for line in fh:
            a, b, s = line.split()
            if int(s) >= MIN_SCORE:
                G.add_edge(a, b, w=int(s))
    print(f"graph: {G.number_of_nodes()} proteins, {G.number_of_edges()} interactions (score >= {MIN_SCORE})")

    # --- Annotate nodes
    func_names = list(FUNCTIONS)
    sys_names = sorted(set(TISSUE_TO_SYSTEM.values())) + [UBIQUITOUS]
    nodes = {}
    for sid in G.nodes:
        name, annot = info.get(sid, (sid.split(".")[1], ""))
        up = uniprot.get(sid) or uniprot_by_gene.get(name)  # UniProt lists one STRING id; fall back to gene
        kws = set(up["Keywords"].split(";")) if up else set()
        acc = up["Entry"] if up else ""
        hpa = hpa_by_acc.get(acc) or hpa_by_gene.get(name)
        classes = set(hpa["Protein class"].split(", ")) if hpa else set()

        funcs = [i for i, f in enumerate(func_names)
                 if kws & FUNCTIONS[f]["kw"] or classes & FUNCTIONS[f]["hpa"]]

        systems, tissues = set(), []
        if hpa:
            level = hpa["RNA tissue specificity"]
            if level in SPECIFIC_LEVELS:
                for t in filter(None, hpa["RNA tissue specific nTPM"].split(";")):
                    tname = t.split(":")[0].strip()
                    tissues.append(tname)
                    if tname in TISSUE_TO_SYSTEM:
                        systems.add(TISSUE_TO_SYSTEM[tname])
            elif level == "Low tissue specificity":
                systems.add(UBIQUITOUS)
            if hpa["RNA blood cell specificity"] in {"Immune cell enriched", "Group enriched"}:
                systems.add("Immune & blood")
        dterms = [d for d in hpa["Disease involvement"].split(", ") if d] if hpa else []
        dset = set(dterms)
        diseases = [i for i, d in enumerate(DISEASES) if dset & DISEASES[d]]
        pnames = parse_protein_names(up["Protein names"]) if up else []
        desc = (pnames[0] if pnames else "") or annot.split(";")[0]
        synonyms = [g.strip() for g in hpa["Gene synonym"].split(",")] if hpa and hpa["Gene synonym"] else []
        aliases = [a for a in dict.fromkeys(pnames[1:] + synonyms) if a and a != name and a != desc]
        nodes[sid] = {
            "name": name,
            "desc": desc,
            "acc": acc,
            "f": funcs,
            "s": sorted(sys_names.index(x) for x in systems),
            "t": tissues,
            "about": annot[:400],
            "ensg": hpa["Ensembl"] if hpa else "",
            "alias": aliases[:12],
            "d": diseases,
            "pf": pfam_of.get(acc, []),
            "dz": [d for d in dterms if d not in HIDDEN_DISEASE_TERMS],
        }

    # --- Communities (modules) for angular placement
    comms = nx.community.louvain_communities(G, weight="w", resolution=1.0, seed=42)
    comms = sorted(comms, key=len, reverse=True)
    small = [c for c in comms if len(c) < 15]
    comms = [c for c in comms if len(c) >= 15] + ([set().union(*small)] if small else [])
    comm_of = {n: i for i, c in enumerate(comms) for n in c}
    print(f"{len(comms)} modules")

    # --- Radial layout: rank by degree -> radius (hubs at centre),
    #     module -> angular sector, neighbours pulled together inside a sector.
    deg = dict(G.degree())
    order = sorted(G.nodes, key=lambda n: -deg[n])
    N = len(order)
    rank = {n: i for i, n in enumerate(order)}
    total = sum(len(c) for c in comms)
    rnd = random.Random(7)
    pos, start = {}, 0.0
    sector = []
    for c in comms:
        span = 2 * math.pi * len(c) / total
        sector.append((start, span))
        start += span
    # initial angle: random inside the sector
    theta = {n: sector[comm_of[n]][0] + rnd.random() * sector[comm_of[n]][1] for n in G.nodes}
    # smooth: a few passes moving each node toward the circular mean of its in-module neighbours
    for _ in range(12):
        new = {}
        for n in G.nodes:
            a0, sp = sector[comm_of[n]]
            nb = [m for m in G[n] if comm_of[m] == comm_of[n]]
            if not nb:
                new[n] = theta[n]
                continue
            x = sum(math.cos(theta[m]) for m in nb) + math.cos(theta[n])
            y = sum(math.sin(theta[m]) for m in nb) + math.sin(theta[n])
            t = math.atan2(y, x) % (2 * math.pi)
            # keep inside the sector, plus jitter so nodes don't collapse
            rel = ((t - a0) % (2 * math.pi))
            rel = min(max(rel, 0.02 * sp), 0.98 * sp)
            new[n] = a0 + 0.6 * rel + 0.4 * (theta[n] - a0) + (rnd.random() - 0.5) * 0.08 * sp
        theta = new
    # The smoothing pulls each module into a narrow spoke. Keep its angular ORDER (neighbours stay close)
    # but spread the nodes evenly over the sector, separately within radius bands so every ring is filled.
    bands = 12
    groups = defaultdict(list)
    for n in G.nodes:
        groups[(comm_of[n], min(bands - 1, int(bands * math.sqrt((rank[n] + 0.5) / N))))].append(n)
    for (c, _), members in groups.items():
        a0, sp = sector[c]
        members.sort(key=lambda n: (theta[n] - a0) % (2 * math.pi))
        for k, n in enumerate(members):
            theta[n] = a0 + sp * (k + 0.5 + (rnd.random() - 0.5) * 0.6) / len(members)
    R = 1000
    for n in G.nodes:
        r = R * math.sqrt((rank[n] + 0.5) / N)  # equal-area rings: top hubs in the middle
        pos[n] = (round(r * math.cos(theta[n]), 1), round(r * math.sin(theta[n]), 1))

    # --- Label each module by its most over-represented function
    fcount_all = Counter(f for n in nodes.values() for f in n["f"])
    modules = []
    for i, c in enumerate(comms):
        fc = Counter(f for n in c for f in nodes[n]["f"])
        best = max(fc, key=lambda f: (fc[f] / len(c)) / (fcount_all[f] / N) * min(1, fc[f] / 5), default=None)
        hub = max(c, key=lambda n: deg[n])
        modules.append({"size": len(c), "hub": nodes[hub]["name"],
                        "label": func_names[best] if best is not None else "Mixed",
                        "a0": sector[i][0], "span": sector[i][1]})

    pf_ids = sorted({x for n in nodes.values() for x in n["pf"]})
    pf_index = {x: i for i, x in enumerate(pf_ids)}

    # --- 3did: Pfam domain pairs seen in contact in PDB structures
    ddi = {}
    with gzip.open(RAW / "3did" / "3did_flat.gz", "rt", errors="replace") as fh:
        cur = None
        for line in fh:
            if line.startswith("#=ID"):
                accs = re.findall(r"(PF\d+)\.\d+@Pfam", line)
                cur = None
                if len(accs) == 2 and accs[0] in pf_index and accs[1] in pf_index:
                    i0, i1 = pf_index[accs[0]], pf_index[accs[1]]
                    swap = i0 > i1  # store the example's chains in (smaller index, larger index) order
                    cur = ddi.setdefault((min(i0, i1), max(i0, i1)), [0, "", "", ""])
            elif line.startswith("#=3D") and cur is not None:
                f = line.split("\t")
                if f[2].split(":")[0] == f[3].split(":")[0]:
                    continue  # same chain: intramolecular domain packing, not a protein-protein interface
                cur[0] += 1
                if not cur[1]:  # example structure: "chain:start-end" for each domain
                    cur[1:] = [f[1], f[3], f[2]] if swap else [f[1], f[2], f[3]]
    print(f"3did domain pairs among network domains: {len(ddi)}")

    # --- Reactome pathways (all levels: a protein in "Insulin receptor signalling" is also in "Signal Transduction")
    rdir = RAW / "reactome"
    pw_name, parent = {}, defaultdict(set)
    with open(rdir / "ReactomePathways.txt", encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[2] == "Homo sapiens":
                pw_name[f[0]] = f[1]
    with open(rdir / "ReactomePathwaysRelation.txt", encoding="utf-8") as fh:
        for line in fh:
            a_, b_ = line.split()
            if a_ in pw_name and b_ in pw_name:
                parent[b_].add(a_)
    acc_to_sid = {nodes[sid]["acc"]: sid for sid in order if nodes[sid]["acc"]}
    members = defaultdict(set)
    with open(rdir / "UniProt2Reactome_All_Levels.txt", encoding="utf-8") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6 and f[5] == "Homo sapiens" and f[1] in pw_name:
                acc = f[0].split("-")[0]
                if acc in acc_to_sid:
                    members[f[1]].add(acc_to_sid[acc])
    pw_keep = sorted((p_ for p_, m in members.items() if len(m) >= 5), key=lambda p_: (-len(members[p_]), pw_name[p_]))
    pw_index = {p_: i for i, p_ in enumerate(pw_keep)}

    def root_of(p_, seen=()):
        ps = [q for q in parent.get(p_, ()) if q not in seen]
        return p_ if not ps else root_of(sorted(ps)[0], seen + (p_,))

    node_pw = defaultdict(list)
    for p_ in pw_keep:
        for sid in members[p_]:
            node_pw[sid].append(pw_index[p_])
    # most specific pathways of a protein: none of its other pathways is a child of them
    children = defaultdict(set)
    for c_, ps in parent.items():
        for q in ps:
            children[q].add(c_)
    node_pw_specific = {}
    for sid, lst in node_pw.items():
        mine = {pw_keep[i] for i in lst}
        node_pw_specific[sid] = sorted((pw_index[p_] for p_ in mine if not (children.get(p_, set()) & mine)),
                                       key=lambda i: len(members[pw_keep[i]]))
    print(f"Reactome: {len(pw_keep)} pathways with >= 5 network proteins; "
          f"{sum(1 for sid in order if node_pw.get(sid))} proteins in at least one")

    # --- curated disease variants (UniProt humsavar, pathogenic / likely pathogenic), written to a lazy-loaded file
    AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
           "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
    net_accs = set(acc_to_sid)
    var_dis, var_dis_idx, variants = [], {}, defaultdict(list)
    with open(RAW / "uniprot_variants" / "humsavar.txt", encoding="utf-8") as fh:
        for line in fh:
            f = line.split()
            if len(f) < 6 or f[1] not in net_accs or f[4] != "LP/P":
                continue
            m = re.fullmatch(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", f[3])
            if not m or m.group(1) not in AA3 or m.group(3) not in AA3:
                continue
            dis = re.sub(r"\s*\[MIM:\d+\]", "", " ".join(f[6:])).strip() or "-"
            if dis not in var_dis_idx:
                var_dis_idx[dis] = len(var_dis); var_dis.append(dis)
            variants[f[1]].append([int(m.group(2)), AA3[m.group(1)] + m.group(2) + AA3[m.group(3)], var_dis_idx[dis]])
    var_out = ROOT / "docs" / "variants.json"
    var_out.write_text(json.dumps({"diseases": var_dis, "variants": variants}, separators=(",", ":")))
    print(f"variants: {sum(map(len, variants.values()))} pathogenic missense in {len(variants)} proteins -> {var_out.name} "
          f"({var_out.stat().st_size / 1e6:.1f} MB)")

    ids = {sid: i for i, sid in enumerate(order)}
    out_nodes = []
    for sid in order:
        n = nodes[sid]
        out_nodes.append([n["name"], n["desc"], n["acc"], deg[sid], *pos[sid], comm_of[sid],
                          n["f"], n["s"], n["t"], n["about"], sid.split(".")[1], n["alias"], n["d"], n["dz"],
                          [pf_index[x] for x in n["pf"]], n["ensg"], sorted(node_pw.get(sid, [])),
                          node_pw_specific.get(sid, [])[:12], len(variants.get(n["acc"], []))])
    # --- signed, directed effects from OmniPath (curated signalling: activation / inhibition)
    # per edge [a, b]: bit0 a activates b, bit1 a inhibits b, bit2 b activates a, bit3 b inhibits a
    sid_of_acc = {nodes[sid]["acc"]: sid for sid in order if nodes[sid]["acc"]}
    effect = defaultdict(int)
    with open(RAW / "omnipath" / "interactions.tsv", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sa, sb = sid_of_acc.get(r["source"]), sid_of_acc.get(r["target"])
            if not sa or not sb or sa == sb or not G.has_edge(sa, sb) or r["is_directed"] not in ("True", "1"):
                continue
            # consensus sign when the sources agree; when they disagree both bits are set (context-dependent)
            T = ("True", "1")
            if r["consensus_stimulation"] in T or r["consensus_inhibition"] in T:
                act, inh = r["consensus_stimulation"] in T, r["consensus_inhibition"] in T
            else:
                act, inh = r["is_stimulation"] in T, r["is_inhibition"] in T
            if not (act or inh):
                continue
            ia, ib = ids[sa], ids[sb]
            key, fwd = (min(ia, ib), max(ia, ib)), ia < ib
            effect[key] |= (act * 1 + inh * 2) << (0 if fwd else 2)
    edges = []
    for a, b, d in G.edges(data=True):
        ia, ib = ids[a], ids[b]
        if ia > ib:
            ia, ib = ib, ia
        m = effect.get((ia, ib), 0)
        edges.append([ia, ib, d["w"], m] if m else [ia, ib, d["w"]])
    print(f"OmniPath: {len(effect)} interactions with a documented activation/inhibition")

    data = {
        "meta": {"source": "STRING v12.0 physical interactions; UniProt; Human Protein Atlas",
                 "min_score": MIN_SCORE, "proteins": N, "interactions": len(edges)},
        "fields": ["name", "desc", "acc", "degree", "x", "y", "module", "functions", "systems",
                   "tissues", "about", "ensp", "aliases", "diseases", "disease_terms", "domains",
                   "ensg", "pathways", "pathways_specific", "n_variants"],
        "functions": func_names,
        "systems": sys_names,
        "diseases": list(DISEASES),
        "domains": [[x, *pfam_info[x]] for x in pf_ids],  # [Pfam id, short name, description]
        # [dom a, dom b, #structures, example PDB, "chain:range" of dom a, "chain:range" of dom b]
        "ddi": [[a, b, n, pdb, ra, rb] for (a, b), (n, pdb, ra, rb) in sorted(ddi.items()) if n],
        # [Reactome id, name, index of its top-level pathway]
        "pathways": [[p_, pw_name[p_], pw_index.get(root_of(p_), -1)] for p_ in pw_keep],
        "modules": modules,
        "nodes": out_nodes,
        "edges": edges,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(data, separators=(",", ":")))
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e6:.1f} MB)")

    # --- coverage report
    fc = Counter(f for n in nodes.values() for f in n["f"])
    sc = Counter(s for n in nodes.values() for s in n["s"])
    print("no function:", sum(1 for n in nodes.values() if not n["f"]),
          "| no system:", sum(1 for n in nodes.values() if not n["s"]))
    for i, f in enumerate(func_names):
        print(f"  F {fc[i]:5d} {f}")
    for i, s in enumerate(sys_names):
        print(f"  S {sc[i]:5d} {s}")
    dc = Counter(d for n in nodes.values() for d in n["d"])
    for i, d in enumerate(DISEASES):
        print(f"  D {dc[i]:5d} {d}")
    used = set().union(*DISEASES.values())
    left = Counter(t for n in nodes.values() for t in n["dz"] if t not in used)
    print("unmapped disease terms:", left.most_common())
    explained = sum(1 for a, b in G.edges if any(
        ddi.get((min(x, y), max(x, y)), [0])[0] > 0
        for x in {pf_index[d] for d in nodes[a]["pf"]} for y in {pf_index[d] for d in nodes[b]["pf"]}))
    print(f"interactions with a known domain-domain contact: {explained} / {G.number_of_edges()}")
    print("proteins with Pfam domains:", sum(1 for n in nodes.values() if n["pf"]), "| families:", len(pf_ids))
    print("top hubs:", [nodes[n]["name"] for n in order[:15]])


if __name__ == "__main__":
    main()
