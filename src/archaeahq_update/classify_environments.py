#!/usr/bin/env python3
"""
Classify archaeal genome assemblies into environmental categories
based on NCBI metadata fields.

Categories (18 total):
 1. Terrestrial_Hot_Spring
 2. Deep_Sea_Hydrothermal_Vent
 3. Cold_Seep
 4. Marine_Water_Column
 5. Marine_Sediment
 6. Coastal_Estuarine
 7. Mangrove
 8. Freshwater
 9. Groundwater_Subsurface_Cave
10. Terrestrial_Soil
11. Peatland_Bog_Wetland
12. Hypersaline_Salt
13. Acid_Mine_Drainage
14. Animal_Host_Associated
15. Wastewater_Biogas_Landfill
16. Oil_Gas_Petroleum
17. Fermentation_Industrial
18. Unknown
"""

import csv
import re
import sys
from collections import Counter

# ---------------------------------------------------------------------------
# Rule definitions
# Each rule is a (category, [positive_patterns], [negative_patterns]) tuple.
# Patterns are case-insensitive substring / regex matches applied to the
# combined metadata string.  The FIRST matching rule wins.
# ---------------------------------------------------------------------------

RULES = [
    # -----------------------------------------------------------------------
    # 1. Terrestrial hot spring  (must be matched BEFORE generic marine/sediment)
    # -----------------------------------------------------------------------
    ("Terrestrial_Hot_Spring", [
        r"\bhot[_\s-]?springs?\b",
        r"\bgeothermal\b",
        r"\bthermophilic\b",
        r"\bboiling\s+springs?\b",
        r"\bthermal\s+springs?\b",
        r"\bthermoacidophilic\b",
        r"\bsolfatara\b",
        r"\bfumarole\b",
        r"\bgeyser\b",
        r"yellowstone",
        r"\bterrestrial.*thermal\b",
        r"\bthermal.*terrestrial\b",
        r"springs?\s+sediment",    # most "spring sediment" is hot spring
        r"terrestrial hot",
        r"\bneutral hot\b",
        r"\balkaline hot\b",
        r"\bacidic hot spring\b",
        r"\bthermophilic sediment\b",
        r"\bhot\s+springs?\s+metagenome\b",
        r"\bmicrobial\s+mat.*terrestrial\b",
        r"\bterrestrial.*microbial\s+mat\b",
    ], [
        r"deep.?sea.*hydrothermal",
        r"hydrothermal\s+vent",
        r"submarine.*hydrothermal",
        r"mid.?ocean.*ridge",
        r"seafloor.*hydrothermal",
        r"\bblack\s+smoker\b",
        r"\bwhite\s+smoker\b",
        r"hydrothermal\s+plume",
    ]),

    # -----------------------------------------------------------------------
    # 2. Deep-sea hydrothermal vent
    # -----------------------------------------------------------------------
    ("Deep_Sea_Hydrothermal_Vent", [
        r"\bhydrothermal\s+vent\b",
        r"\bhydrothermal\s+chimney\b",
        r"\bblack\s+smoker\b",
        r"\bwhite\s+smoker\b",
        r"\bhydrothermal\s+deposit\b",
        r"deep.?sea.*hydrothermal",
        r"marine\s+hydrothermal",
        r"hydrothermal\s+plume",
        r"mid.?ocean.*ridge",
        r"mid.?atlantic.*ridge",
        r"east\s+pacific\s+rise",
        r"seafloor.*hydrothermal",
        r"submarine\s+hydrothermal",
        r"hydrothermal\s+field",
        r"hydrothermal\s+system.*marine",
        r"vent\s+fluid",
        r"\bhydrothermal\s+metagenome\b",
        r"hydrothermal.*sediment.*marine",
        r"ocean.*vent",
    ], []),

    # -----------------------------------------------------------------------
    # 3. Cold seep / methane seep
    # -----------------------------------------------------------------------
    ("Cold_Seep", [
        r"\bcold\s+seep\b",
        r"\bmethane\s+seep\b",
        r"\bgas\s+seep\b",
        r"\bgas\s+hydrate\b",
        r"\bhydrate\s+seep\b",
        r"\bmethane\s+hydrate\b",
        r"\bseep\s+sediment\b",
        r"\bseep\b.*\bsediment\b",
        r"\bcold\s+seep\s+metagenome\b",
        r"\basphalt\s+volcano\b",
        r"\bmud\s+volcano\b",
    ], []),

    # -----------------------------------------------------------------------
    # 4. Marine water column / open ocean
    # -----------------------------------------------------------------------
    ("Marine_Water_Column", [
        r"\bmarine\s+metagenome\b",
        r"\bseawater\s+metagenome\b",
        r"\bmarine\s+plankton\b",
        r"\bpelagic\b",
        r"\bocean\s+water\b",
        r"\bopen\s+ocean\b",
        r"\bdeep\s+ocean\b",
        r"\bocean\s+metagenome\b",
        r"\bmarine\s+water\b",
        r"\bdeep\s+atlantic\b",
        r"\bdeep\s+pacific\b",
        r"\btara\s+ocean",
        r"\bgo.ship\b",
        r"\bwater\s+column\b",
        r"\boxygen\s+minimum\s+zone\b",
        r"\bhypoxic\s+seawater\b",
        r"\banoxic\s+seawater\b",
        r"\bsea\s+water\b",
        r"\b(sea|ocean|marine)\s+surface\s+water\b",
        r"\bmarine\s+water\s+sample\b",
        r"\bcoastal\s+(sea\s*)?water\b",
        r"\bbay\s+water\b",
        r"\bphotic\s+zone\b",
        r"north sea",
        r"atlantic ocean.*water",
        r"pacific ocean.*water",
        r"indian ocean.*water",
        r"arctic ocean.*water",
        r"\bnorth\s+atlantic\b",
        r"\bsouth\s+atlantic\b",
        r"\barctic.*sea\b",
        r"\bsea.?ice\b",
        r"\bocean.*epipelagic\b",
        r"\bross\s+ice\s+shelf\b",
        r"beneath.*ice\s+shelf",
        r"\bseawater\b",
    ], [
        r"\bsediment\b(?!\s+trap)",     # exclude sediment but not "sediment trap"
        r"\bmangrove\b",
        r"\bcoastal.*sediment\b",
        r"\bestuar",
        r"salt marsh",
        r"\btidal\b",
        r"intertidal",
    ]),

    # -----------------------------------------------------------------------
    # 5. Marine sediment (deep-sea and coastal marine, excl. hydrothermal/seep)
    # -----------------------------------------------------------------------
    ("Marine_Sediment", [
        r"\bmarine\s+sediment\b",
        r"\bdeep.?sea\s+sediment\b",
        r"\bdeep\s+sea\s+sediment\b",
        r"\bmarine\s+sediment\s+metagenome\b",
        r"\bbenthic\b",
        r"\babyssal\b",
        r"\bseafloor\s+sediment\b",
        r"\bocean\s+(floor|bed|bottom)\b",
        r"\bdeep\s+sea\s+metagenome\b",
        r"\bbottom\s+of.*ocean\b",
        r"\bsea\s*floor\b",
        r"\bsea\s+floor\b",
        r"\bsubseafloor\b",
        r"\bsub.?seafloor\b",
    ], [
        r"\bhydrothermal\b",
        r"\bcold\s+seep\b",
        r"\bmethane\s+seep\b",
        r"\bmangrove\b",
        r"\bestuar",
        r"\bintertidal\b",
        r"\bsalt\s+marsh\b",
        r"\btidal\s+flat\b",
    ]),

    # -----------------------------------------------------------------------
    # 6. Coastal / Estuarine
    # -----------------------------------------------------------------------
    ("Coastal_Estuarine", [
        r"\bestuar",
        r"\bintertidal\b",
        r"\btidal\s+flat\b",
        r"\btidal\s+sand",
        r"\bmudflat\b",
        r"\bsalt\s+marsh\b",
        r"\bcoastal\s+sediment\b",
        r"\bcoastal\s+marine\b",
        r"\bcoastal\s+bare\b",
        r"\bnear.?shore\b",
        r"\bcoastal\s+metagenome\b",
        r"estuary metagenome",
        r"\bbeach\b",
        r"\blagoon\b",
        r"\btidal\s+zone\b",
        r"\bshore\b",
        r"riparian\s+sediment",
        r"\bshoreline\b",
        r"\bcoastal\b",
        r"\bfloodplain\b",
        r"\bcoastal.*seawater\b",
        r"\bdelta\b",
        r"\bblack\s+sea\b",
        r"\bcaspian\s+sea\b",
    ], [
        r"\bmangrove\b",
        r"\bsalt\s+lake\b",
        r"\bhypersaline\b",
        r"\bsalt\s+mine\b",
    ]),

    # -----------------------------------------------------------------------
    # 7. Mangrove
    # -----------------------------------------------------------------------
    ("Mangrove", [
        r"\bmangrove\b",
    ], []),

    # -----------------------------------------------------------------------
    # 8. Freshwater (lakes, rivers, streams, ponds — non-saline)
    # -----------------------------------------------------------------------
    ("Freshwater", [
        r"\bfreshwater\b",
        r"\blake\s+water\b",
        r"\blake\s+sediment\b",
        r"\blake\s+metagenome\b",
        r"\bfreshwater\s+sediment\b",
        r"\bfreshwater\s+metagenome\b",
        r"\bfreshwater\s+lake\b",
        r"\briver\s+water\b",
        r"\briver\s+sediment\b",
        r"\bstream\b",
        r"\bfreshwater\s+pond\b",
        r"\bpond\s+sediment\b",
        r"\blentic\b",
        r"\blotic\b",
        r"\bmeromictic\s+lake\b",
        r"\bfrozen\s+lake\b",
        r"\bperiphyton\b",
        r"\blake\s+water\b",
        r"\bwater\s+reservoir\b",
        r"freshwater.*microbial mat",
    ], [
        r"\bsaline\b",
        r"\bhypersaline\b",
        r"\bsalt\s+lake\b",
        r"\bsoda\s+lake\b",
        r"\balkaline\s+lake\b",
        r"\bhot\s+spring\b",
    ]),

    # -----------------------------------------------------------------------
    # 9. Groundwater / Aquifer / Cave / Deep subsurface
    # -----------------------------------------------------------------------
    ("Groundwater_Subsurface_Cave", [
        r"\bgroundwater\b",
        r"\bground\s+water\b",
        r"\baquifer\b",
        r"\bsubsurface\b",
        r"\bborehole\b",
        r"\bunderground\s+water\b",
        r"\bdeep\s+subsurface\b",
        r"\bdeep\s+aquifer\b",
        r"\bsubterranean\b",
        r"aspo\s+hrl",
        r"\bolkiluoto\b",
        r"\bcave\b",
        r"\bkarst\b",
        r"\bsulfidic\s+cave\b",
        r"\bsulfidic\s+groundwater\b",
        r"\bmine\s+water\b",
        r"\bdeep\s+rock\b",
        r"\brock\s+borehole\b",
        r"\bdeep\s+coal\s+seam\b",
        r"\bcoal\s+seam\b",
        r"\brock\s+metagenome\b",
        r"\bdeep\s+rock\s+metagenome\b",
        r"\bbasaltic\s+crust\b",
        r"\bbasaltic\s+aquifer\b",
        r"\bbasalt.*fluid\b",
        r"\bcrustal\s+fluid\b",
        r"\bbracki(sh)?\s+groundwater\b",
    ], [
        r"\bmine\s+drain",
        r"\bacid\s+mine\b",
        r"\bAMD\b",
        r"\boil\s+(field|well|reservoir)\b",
        r"\bgas\s+(well|field)\b",
    ]),

    # -----------------------------------------------------------------------
    # 10. Terrestrial soil (incl. permafrost, agricultural, desert, forest, arctic)
    # -----------------------------------------------------------------------
    ("Terrestrial_Soil", [
        r"\bsoil\b",
        r"\bsoil\s+metagenome\b",
        r"\bpaddy\s+field\b",
        r"\bpaddy\s+soil\b",
        r"\bcropland\b",
        r"\bagricultural\s+soil\b",
        r"\bfarm\s+soil\b",
        r"\bgrassland\s+soil\b",
        r"\bgrassland\b",
        r"\bforest\s+soil\b",
        r"\bpermafrost\b",
        r"\bactive\s+layer\s+soil\b",
        r"\btundra\b",
        r"\bdesert\s+soil\b",
        r"\bdesert\s+biocrust\b",
        r"\barid\s+soil\b",
        r"\brhizosphere\b",
        r"\bcompost\b",
        r"\btopsoil\b",
        r"\bterrestrial\s+metagenome\b",
        r"\bsaline\s+soil\b",
        r"\bsaltern\s+soil\b",
        r"\bsalted\s+soil\b",
        r"\bsalt\s+flat\b",       # salar / playa
        r"\bsalar\b",
        r"\bsaline\s+lake.*soil\b",
        r"\bsediment\s+metagenome.*terrestrial\b",
        r"\bterrestrial.*sediment\s+metagenome\b",
        r"\bmetagenome.*terrestrial\b",
        r"\bterrestrial.*metagenome\b",
        r"\bsteppe\b",
        r"\bprairie\b",
        r"\bshrubland\b",
        r"\bterrestrial\s+biome\b",
        r"\bantarct.*soil\b",
        r"\balpine\s+soil\b",
        r"\bvolcanic\s+soil\b",
    ], [
        r"\bhot\s+spring\b",
        r"\bpeat\b",
        r"\bbog\b",
        r"\bwetland\b",
        r"\bmangrove\b",
    ]),

    # -----------------------------------------------------------------------
    # 11. Peatland / Bog / Wetland
    # -----------------------------------------------------------------------
    ("Peatland_Bog_Wetland", [
        r"\bpeat\b",
        r"\bbog\b",
        r"\bfen\b",
        r"\bmire\b",
        r"\bwetland\b",
        r"\bmarsh\b",
        r"\bswamp\b",
        r"\bsphagnum\b",
        r"\bpeatland\b",
        r"\bpalsa\b",
    ], [
        r"\bsalt\s+marsh\b",
    ]),

    # -----------------------------------------------------------------------
    # 12. Hypersaline / Salt environments
    # -----------------------------------------------------------------------
    ("Hypersaline_Salt", [
        r"\bhypersaline\b",
        r"\bsalt\s+lake\b",
        r"\bsalt\s+mine\b",
        r"\bsolar\s+salt\b",
        r"\bsalt\s+crust\b",
        r"\bhalite\b",
        r"\bsoda\s+lake\b",
        r"\bsaline\s+lake\b",
        r"\bsaline\s+evaporation\s+pond\b",
        r"\bsaltern\b",
        r"\bevaporation\s+pond\b",
        r"\bsalt\s+pond\b",
        r"\bhalo(phile|philic|archaea)\b",
        r"\bhypersaline\s+lake\b",
        r"\bsalt\s+pan\b",
        r"\bdead\s+sea\b",
        r"\bgreat\s+salt\s+lake\b",
        r"\bhighly\s+saline\b",
        r"\bsalt\s+flat\b",
        r"\bsalt.*brine\b",
        r"\bsaline\s+water\b",
        r"\bbrine\b",
        r"\brock\s+salt\b",
        r"\bsubterranean.*salt\b",
        r"commercial salt",
        r"solar salt",
        r"\bsalt\s+work",
        r"\bhalophile\b",
        r"\bedible\s+salt\b",
        r"\b(lake|laguna)\s+(lejia|magadi|asal|assal|natron|bogoria)\b",
        r"\bsaline\s+sediment\b",
        r"\bsalt\s+lake\s+sediment\b",
        r"\bendorheic\b",
        r"(?<!\w)salt(?!\s+marsh)(?!\s+water)(?!\s+flat)(?!\s*water)",  # "salt" alone

    ], [
        r"\bsalt\s+marsh\b",
        r"\bmine\s+drain",
        r"\bacid\s+mine\b",
    ]),

    # -----------------------------------------------------------------------
    # 13. Acid mine drainage / Mine
    # -----------------------------------------------------------------------
    ("Acid_Mine_Drainage", [
        r"\bacid\s+mine\s+drain",
        r"\bamd\b",
        r"\bmine\s+drain",
        r"\bmine\s+tailings\b",
        r"\btailing(s)?\s+pond\b",
        r"\bmine\s+water\b",
        r"\brichmond\s+mine\b",
        r"\bcopper\s+mine\b",
        r"\byellow\s+boy\b",
        r"\bbioleaching\b",
        r"\bmine\s+drainage\s+metagenome\b",
        r"\bmineral\s+leaching\b",
        r"\bacid\s+rock\s+drainage\b",
        r"\bpyrite\s+sediment\b",
        r"\bfankou\b",
        r"mine.*acidic",
        r"acidic.*mine",
        r"\bsulfuric\s+acid.*mine\b",
        r"\bpyrite\s+surface\b",
        r"\bslime\s+streamer\b",
        r"\bsulfide\s+ore\b",
    ], []),

    # -----------------------------------------------------------------------
    # 14. Animal / Host-associated (gut, rumen, feces, body fluids, sponge, coral)
    # -----------------------------------------------------------------------
    ("Animal_Host_Associated", [
        r"\brumen\b",
        r"\bgut\b",
        r"\bfece(s)?\b",
        r"\bfaece(s)?\b",
        r"\bfecal\b",
        r"\bfaecal\b",
        r"\bgastrointestinal\b",
        r"\bintestinal\b",
        r"\bintestine\b",
        r"\bcaecum\b",
        r"\bcecum\b",
        r"\bhuman\s+gut\b",
        r"\bbovine\s+gut\b",
        r"\bpig\s+gut\b",
        r"\btermite\s+gut\b",
        r"\bmarine\s+sponge\b",
        r"\bsponge\s+(tissue|metagenome)\b",
        r"\bcoral\b",
        r"\bholothurian\b",
        r"\bsea\s+cucumber\b",
        r"\bhost.?associated\b",
        r"\bbodily\s+fluid\b",
        r"\banimal\s+gut\b",
        r"\banimal.*feces\b",
        r"\bruminant\b",
        r"\bcattle\s+rumen\b",
        r"\bfecal\s+sample\b",
        r"\bfecal\s+metagenome\b",
        r"\bfeces\s+metagenome\b",
        r"\bgut\s+metagenome\b",
        r"\bciliate\b",
        r"\bgastrointestinal\b",
        r"\bcaecal\b",
        r"\bcecal\b",
        r"gastrointestinal_tract",
        r"\bcaecal_contents\b",
        r"\bannelid\b",
        r"\bannelid\s+tissue\b",
    ], []),

    # -----------------------------------------------------------------------
    # 15. Wastewater / Biogas / Landfill / Bioreactor (engineered systems)
    # -----------------------------------------------------------------------
    ("Wastewater_Biogas_Landfill", [
        r"\bwaste\s*water\b",
        r"\bwastewater\b",
        r"\bactivated\s+sludge\b",
        r"\bbiogas\b",
        r"\banaerobic\s+dige[rs]t",
        r"\banaerobic\s+bioreactor\b",
        r"\bbioreactor\b",
        r"\bsewage\b",
        r"\bsludge\b",
        r"\blandfill\b",
        r"\bleachate\b",
        r"\bpulp.*mill\b",
        r"\bpaper.*mill\b",
        r"\bpalm\s+oil\s+mill\b",
        r"\bWWTP\b",
        r"\bMBBR\b",
        r"\bAnMBR\b",
        r"\baeration\s+tank\b",
        r"\bnitrification\b",
        r"\bnitritation\b",
        r"\banammox\b",
        r"\bbiofilm.*reactor\b",
        r"\banaerobic\s+reactor\b",
        r"\bwaste\s+treatment\b",
        r"\bindustrial\s+waste\s+metagenome\b",
        r"\bengineered\b",
    ], [
        r"\bhot\s+spring\b",
        r"\bhydrothermal\b",
    ]),

    # -----------------------------------------------------------------------
    # 16. Oil / Gas / Petroleum subsurface
    # -----------------------------------------------------------------------
    ("Oil_Gas_Petroleum", [
        r"\boil\s+field\b",
        r"\boil\s+well\b",
        r"\boil\s+reservoir\b",
        r"\bgas\s+field\b",
        r"\bgas\s+well\b",
        r"\bhydrocarbon\b",
        r"\bpetroleum\b",
        r"\bproduced\s+(water|fluid)\b",
        r"\bformation\s+water\b",
        r"\bhydraulic\s+fract",
        r"\bshale\s+gas\b",
        r"\bfracking\b",
        r"\boil\s+sands\b",
        r"\btar\s+sand\b",
        r"\boil\s+shale\b",
        r"\bpetroleum\s+reservoir\b",
        r"\boil\s+production\b",
        r"\bgas\s+reservoir\b",
        r"\bnatural\s+gas\s+well\b",
        r"\bcorroded\s+pipe\b",
        r"oil\s+(field|production)\s+metagenome",
        r"\bhydrocarbon\s+metagenome\b",
        r"\bnaphtha\b",
        r"\bbitumen\b",
        r"\brifle.*sediment\b",
        r"\brifle.*well\b",
    ], []),

    # -----------------------------------------------------------------------
    # 17. Fermentation / Industrial food / Aquaculture
    # -----------------------------------------------------------------------
    ("Fermentation_Industrial", [
        r"\bferment",
        r"\bbaijiu\b",
        r"\bkimchi\b",
        r"\bmiso\b",
        r"\bfood\s+product\b",
        r"\baquaculture\b",
        r"\bfishpond\b",
        r"\bfish\s+pond\b",
        r"\bbioelectrochemical\b",
        r"\bbiofuel\b",
        r"\bballast\s+tank\b",
        r"\bferment.*metagenome\b",
        r"\bbiooxidation\b",
        r"\bfood\s+production\b",
        r"\bfood\s+product\s+metagenome\b",
    ], [
        r"\banaerobic\s+dige[rs]t",
        r"\bbiogas\b",
    ]),
]

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def normalise(text: str) -> str:
    """Lower-case, collapse whitespace."""
    return re.sub(r"\s+", " ", str(text).lower().strip())


def classify(row: dict) -> str:
    """Return the first matching category label, else 'Unknown'."""
    # Build combined text from all metadata fields (excluding Assembly ID)
    parts = [
        row.get("Metagenome_soource", ""),
        row.get("geo_loc_name", ""),
        row.get("Isolation_Source", ""),
        row.get("Env_Broad", ""),
        row.get("Env_Local", ""),
        row.get("Env_medium", ""),
    ]
    combined = normalise(" | ".join(str(p) for p in parts))

    # Skip completely empty / missing records immediately
    empty_tokens = {"0", "missing", "n/a", "not applicable", "not collected",
                    "not available", "na", "none", ""}
    clean_parts = [p.strip().lower() for p in parts if p.strip().lower() not in empty_tokens]
    if not clean_parts:
        return "Unknown"

    for category, pos_patterns, neg_patterns in RULES:
        # Check positive patterns
        matched = any(re.search(p, combined, re.IGNORECASE) for p in pos_patterns)
        if not matched:
            continue
        # Check negative (exclusion) patterns
        excluded = any(re.search(p, combined, re.IGNORECASE) for p in neg_patterns)
        if excluded:
            continue
        return category

    return "Unknown"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    input_file  = "Info_unify-Sample_environment.txt"
    output_file = "Classified_environments.tsv"

    with open(input_file, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)

    counts: Counter = Counter()
    output_rows = []
    for row in rows:
        label = classify(row)
        counts[label] += 1
        output_rows.append({**row, "Env_Category": label})

    # Write classified output
    fieldnames = list(rows[0].keys()) + ["Env_Category"]
    with open(output_file, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    # Summary table
    summary_file = "Env_Category_Summary.tsv"
    total = len(rows)
    summary_rows = sorted(counts.items(), key=lambda x: -x[1])

    with open(summary_file, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["Env_Category", "Count", "Percentage"])
        for cat, n in summary_rows:
            writer.writerow([cat, n, f"{100 * n / total:.2f}"])
        writer.writerow(["TOTAL", total, "100.00"])

    # Print to stdout
    print(f"\n{'Category':<35} {'Count':>6}  {'%':>6}")
    print("-" * 52)
    for cat, n in summary_rows:
        print(f"{cat:<35} {n:>6}  {100*n/total:>5.1f}%")
    print("-" * 52)
    print(f"{'TOTAL':<35} {total:>6}")
    print(f"\nClassified output written to : {output_file}")
    print(f"Category summary written to  : {summary_file}")


if __name__ == "__main__":
    main()
