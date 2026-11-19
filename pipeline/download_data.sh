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
mkdir -p reactome uniprot_variants
for f in UniProt2Reactome_All_Levels.txt ReactomePathways.txt ReactomePathwaysRelation.txt; do
  curl -fSL -o "reactome/$f" "https://reactome.org/download/current/$f"
done
curl -fSL -o uniprot_variants/humsavar.txt https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/humsavar.txt
mkdir -p omnipath
curl -fSL -o omnipath/interactions.tsv "https://omnipathdb.org/interactions?genesymbols=yes&fields=sources,references&license=academic"
curl -fSL -o string_physical_detailed.txt.gz https://stringdb-downloads.org/download/protein.physical.links.detailed.v12.0/9606.protein.physical.links.detailed.v12.0.txt.gz
mkdir -p biogrid && curl -fSL -o biogrid/all.zip https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-ORGANISM-LATEST.tab3.zip \
  && unzip -o -q biogrid/all.zip "*Homo_sapiens*" -d biogrid && rm biogrid/all.zip
echo "Done. Now run: python3 pipeline/build_data.py"
