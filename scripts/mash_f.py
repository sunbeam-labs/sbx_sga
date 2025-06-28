import re
import os


def open_report(report):
    sample_name = os.path.basename(report.split("_sorted_winning.tab")[0])
    with open(report, "r") as report_obj:
        filelines = report_obj.readlines()
    top_lines = filelines[:20]
    return sample_name, top_lines


def process_mash_line(line):
    line_list = line.rstrip().split("\t")
    species_line = line_list[-1]
    matches = re.findall("N[A-Z]_[0-9A-Z]+\.[0-9]", species_line)
    if not matches:
        try:
            matches = re.findall("[A-Z]{2}_[0-9]+\\.[0-9]", species_line)
        except:
            raise ValueError(f"No match found in species_line: {species_line}")
    split_char = matches[0]
    species_split = line.split(split_char)[1].lstrip()
    species = " ".join(species_split.split()[:2])
    median_multiplicity = float(line_list[2])
    identity = float(line_list[0])
    hits = int(line_list[1].split("/")[0])
    return species, median_multiplicity, identity, hits

def get_first_non_phage_hit(lines):
    for idx, line in enumerate(lines):
        if "phage" not in line.lower():
            return process_mash_line(line), idx
    raise ValueError("All top hits contain 'phage' — no valid non-phage hit found.")

def parse_report(top_lines):
    target_species = []

    # Get top non-phage hit and its index
    (top_species, top_median_multiplicity, top_identity, top_hits), top_index = get_first_non_phage_hit(top_lines)   
    
    if (top_identity >= 0.85) and (top_hits >= 100):
        target_species.append(top_species)

    # Set the threshold for median multiplicity
    threshold = 0.05 * top_median_multiplicity

    # Iterate through the rest of the hits, excluding top_index
    for i, line in enumerate(top_lines):
        if i == top_index:
            continue
        species, median_multiplicity, identity, hits = process_mash_line(line)
        if (identity >= 0.85) and (hits >= 100):
            if any(term in species for term in ["phage", "Phage", "sp."]):
                continue
            if median_multiplicity >= threshold:
                target_species.append(species)

    return set(target_species)

def contamination_call(target_set):
    mash_dict = {}
    if len(target_set) <= 1:
        mash_dict["NA"] = ""
    else:
        species = " ".join(sorted(target_set))
        mash_dict["Contaminated"] = species
    return mash_dict


def write_report(output, sample_name, mash_dict):
    # Expecting that the dictionary is just one key-value pair, so need to check that
    if len(mash_dict) == 1:
        status = list(mash_dict.keys())[0]
    else:
        # Raise error if dictionary is not proper length
        raise ValueError(
            f"Expected mash_dict to have exactly one key-value pair, but got {len(mash_dict)}."
        )
    with open(output, "w") as out:
        if status == "Contaminated":
            contaminated_spp = mash_dict[status]
            out.write(f"{sample_name}\tContaminated\t{contaminated_spp}\n")
        else:
            out.write(f"{sample_name}\tNA\tNA\n")
    return output
