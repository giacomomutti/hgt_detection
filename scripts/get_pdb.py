import argparse
from unipressed import IdMappingClient, UniprotkbClient
from time import sleep
import os

def parse_args():
    parser=argparse.ArgumentParser(description="Get a new taxdump including new close proteomes")
    parser.add_argument('-i', '--input', dest='input', required=True, type=str,
                    help='Input code'),
    parser.add_argument('-o', '--outdir', dest='outdir', required=True, type=str,
            help='output folder'),
    args=parser.parse_args()
    return args


if __name__ == "__main__":

    inputs=parse_args()
    
    with open(inputs.input) as infile:
        proteins = set(el.strip() for el in infile)

    if not os.path.exists(inputs.outdir):
        os.makedirs(inputs.outdir)


    request = IdMappingClient.submit(
        source="EMBL-GenBank-DDBJ_CDS", dest="UniProtKB", ids=proteins
    )

    while True:
        status = request.get_status()
        if status in {"FINISHED", "ERROR"}:
            break
        else:
            sleep(3)


    uniprot_ids = {el["to"]:el["from"] for el in list(request.each_result())}

    for uniprot,protein in uniprot_ids.items():
        model_url = f'https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb'
        os.system(f'curl {model_url} -o {inputs.outdir}/{protein}.pdb')
