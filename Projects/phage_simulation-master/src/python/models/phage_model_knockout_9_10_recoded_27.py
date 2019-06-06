from Bio import Entrez, SeqIO
import pinetree as pt

CELL_VOLUME = 1.1e-15
PHI10_BIND = 1.82e7  # Binding constant for phi10

IGNORE_REGULATORY = ["E. coli promoter E[6]",
                     "T7 promoter phiOR",
                     "T7 promoter phiOL",
                     "E. coli promoter A0 (leftward)"]

IGNORE_GENES = ["gene 10B",
                "possible gene 5.5-5.7",
                "gene 1.5",
                "gene 1.6",
                "gene 4.1",
                "gene 4B",
                "gene 0.6A",
                "gene 0.6B",
                "possible gene 0.6B",
                "gene 0.5",
                "gene 0.4"]

RELABEL_GENES = {"gene 2": "gp-2",
                 "gene 1": "rnapol-1",
                 "gene 3.5": "lysozyme-3.5",
                 "gene 0.7": "protein_kinase-0.7"}


def get_promoter_interactions(name):
    '''
    Calculate promoter binding strengths. The relative strengths defined here
    come from 2012 Covert, et al paper.
    '''
    ecoli_strong = ["E. coli promoter A1",
                    "E. coli promoter A2",
                    "E. coli promoter A3"]
    ecoli_weak = ["E. coli B promoter",
                  "E. coli C promoter"]
    phi1_3 = ["T7 promoter phi1.1A",
              "T7 promoter phi1.1B",
              "T7 promoter phi1.3",
              "T7 promoter phi1.5",
              "T7 promoter phi1.6"]
    phi3_8 = ["T7 promoter phi2.5",
              "T7 promoter phi3.8",
              "T7 promoter phi4c",
              "T7 promoter phi4.3",
              "T7 promoter phi4.7"]
    phi6_5 = ["T7 promoter phi6.5"]
    phi9 = ["T7 promoter phi9"]
    phi10 = ["T7 promoter phi10"]
    phi13 = ["T7 promoter phi13",
             "T7 promoter phi17"]

    if name in ecoli_strong:
        return {'ecolipol': 10e4,
                'ecolipol-p': 3e4}
    elif name in ecoli_weak:
        return {'ecolipol': 1e4,
                'ecolipol-p': 0.3e4}
    elif name in phi1_3:
        return {'rnapol-1': PHI10_BIND * 0.01,
                'rnapol-3.5': PHI10_BIND * 0.01 * 0.5}
    elif name in phi3_8:
        return {'rnapol-1': PHI10_BIND * 0.01,
                'rnapol-3.5': PHI10_BIND * 0.01 * 0.5}
    elif name in phi6_5:
        return {'rnapol-1': PHI10_BIND * 0.05,
                'rnapol-3.5': PHI10_BIND * 0.05}
    elif name in phi9:
        return {'rnapol-1': 0,
                'rnapol-3.5': 0}
    elif name in phi10:
        return {'rnapol-1': 0,
                'rnapol-3.5': 0}
    elif name in phi13:
        return {'rnapol-1': PHI10_BIND * 0.1,
                'rnapol-3.5': PHI10_BIND * 0.1}
    else:
        raise ValueError(
            "Promoter strength for {0} not assigned.".format(name))


def get_terminator_interactions(name):
    '''
    Get terminator efficiencies.
    '''
    if name == "E. coli transcription terminator TE":
        return {'ecolipol': 1.0,
                'ecolipol-p': 1.0,
                'rnapol-1': 0.0,
                'rnapol-3.5': 0.0}
    elif name == "T7 transcription terminator Tphi":
        return {'rnapol-1': 0.85,
                'rnapol-3.5': 0.85}
    else:
        return {'name': 0.0}


def main():
    sim = pt.Model(cell_volume=CELL_VOLUME)

    # Download T7 wild-type genbank records
    Entrez.email = "benjamin.r.jack@gmail.com"
    handle = Entrez.efetch(db="nuccore",
                           id=["NC_001604"],
                           rettype="gb",
                           retmode="text")

    record = SeqIO.read(handle, "genbank")
    genome_length = len(record.seq)
    phage = pt.Genome(name="phage", length=genome_length)

    for feature in record.features:
        weights = [1.0] * len(record.seq)
        # Convert to inclusive genomic coordinates
        start = feature.location.start.position + 1
        stop = feature.location.end.position
        name = ''
        if "note" in feature.qualifiers:
            name = feature.qualifiers["note"][0]
        # Grab promoters and terminators
        if feature.type == "regulatory":
            if name in IGNORE_REGULATORY:
                continue
            # Construct promoter
            if "promoter" in feature.qualifiers["regulatory_class"]:
                length = stop - start
                if length < 35:
                    start = start - 35
                interactions = get_promoter_interactions(name)
                phage.add_promoter(name, start, stop, interactions)
            # Construct terminator params
            if "terminator" in feature.qualifiers["regulatory_class"]:
                interactions = get_terminator_interactions(name)
                phage.add_terminator(name, start, stop, interactions)
        # Grab genes/CDSes
        if feature.type == "gene":
            if name in IGNORE_GENES:
                continue
            if name in RELABEL_GENES:
                name = RELABEL_GENES[name]
            # Construct CDS parameters for this gene
            # print(name)
            phage.add_gene(name=name, start=start, stop=stop,
                           rbs_start=start - 30, rbs_stop=start, rbs_strength=1e7)
        # Recode gene 10A
        if name == "gene 10A":
            gene10_start = start
            gene10_stop = stop

    weights[gene10_start:gene10_stop] = [0.27] * (gene10_stop - gene10_start)

    mask_interactions = ["rnapol-1", "rnapol-3.5",
                         "ecolipol", "ecolipol-p", "ecolipol-2", "ecolipol-2-p"]
    phage.add_mask(500, mask_interactions)

    phage.add_weights(weights)

    sim.register_genome(phage)

    sim.add_polymerase("rnapol-1", 35, 230, 0)
    sim.add_polymerase("rnapol-3.5", 35, 230, 0)
    sim.add_polymerase("ecolipol", 35, 45, 0)
    sim.add_polymerase("ecolipol-p", 35, 45, 0)
    sim.add_polymerase("ecolipol-2", 35, 45, 0)
    sim.add_polymerase("ecolipol-2-p", 35, 45, 0)

    sim.add_ribosome(30, 30, 0)

    sim.add_species("bound_ribosome", 10000)

    sim.add_species("bound_ecolipol", 1800)
    sim.add_species("bound_ecolipol_p", 0)
    sim.add_species("ecoli_genome", 0)
    sim.add_species("ecoli_transcript", 0)

    sim.add_reaction(1e6, ["ecoli_transcript", "__ribosome"], [
                     "bound_ribosome"])

    sim.add_reaction(0.04, ["bound_ribosome"], [
                     "__ribosome", "ecoli_transcript"])

    sim.add_reaction(0.001925, ["ecoli_transcript"], ["degraded_transcript"])

    sim.add_reaction(1e7, ["ecolipol", "ecoli_genome"], ["bound_ecolipol"])

    sim.add_reaction(
        0.3e7, ["ecolipol-p", "ecoli_genome"], ["bound_ecolipol_p"])

    sim.add_reaction(0.04, ["bound_ecolipol"], [
                     "ecolipol", "ecoli_genome", "ecoli_transcript"])

    sim.add_reaction(0.04, ["bound_ecolipol_p"], [
                     "ecolipol-p", "ecoli_genome", "ecoli_transcript"])

    sim.add_reaction(3.8e7, ["protein_kinase-0.7", "ecolipol"],
                     ["ecolipol-p", "protein_kinase-0.7"])

    sim.add_reaction(3.8e7, ["protein_kinase-0.7", "ecolipol-2"],
                     ["ecolipol-2-p", "protein_kinase-0.7"])

    sim.add_reaction(3.8e7, ["gp-2", "ecolipol"], ["ecolipol-2"])

    sim.add_reaction(3.8e7, ["gp-2", "ecolipol-p"], ["ecolipol-2-p"])

    sim.add_reaction(1.1, ["ecolipol-2-p"], ["gp-2", "ecolipol-p"])

    sim.add_reaction(1.1, ["ecolipol-2"], ["gp-2", "ecolipol"])

    sim.add_reaction(3.8e9, ["lysozyme-3.5", "rnapol-1"], ["rnapol-3.5"])

    sim.add_reaction(3.5, ["rnapol-3.5"], ["lysozyme-3.5", "rnapol-1"])

    sim.seed(54)

    sim.simulate(time_limit=1200, time_step=5,
                 output="phage_model_knockout_9_10_recoded_27.tsv")


if __name__ == "__main__":
    main()
