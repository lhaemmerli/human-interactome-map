#!/usr/bin/env bash
# Download the raw inputs used by build_data.py into data/raw/.
set -euo pipefail
cd "$(dirname "$0")/../data" && mkdir -p raw/3did && cd raw
UP="https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=tsv&query=organism_id:9606+AND+reviewed:true"
curl -fSL -o string_physical.txt.gz https://stringdb-downloads.org/download/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz
curl -fSL -o string_info.txt.gz      https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz
curl -fSL -o proteinatlas.tsv.zip    https://www.proteinatlas.org/download/proteinatlas.tsv.zip && unzip -o -q proteinatlas.tsv.zip
curl -fSL -o uniprot_human.tsv.gz    "$UP&fields=accession,gene_primary,protein_name,keyword,xref_string,cc_function"
curl -fSL -o uniprot_pfam.tsv.gz     "$UP&fields=accession,xref_pfam"
curl -fSL -o pfam_clans.tsv.gz       https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz
curl -fSL --retry 3 -o 3did/3did_flat.gz https://3did.irbbarcelona.org/download/current/3did_flat.gz
echo "Done. Now run: python3 pipeline/build_data.py"
