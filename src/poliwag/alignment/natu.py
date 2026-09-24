from antismash.modules.nrps_pks.name_mappings import get_substrate_by_name

def convert_to_natu_format(polymer: list[str]) -> list[str]:
    natu_formatted_polymer = []
    for monomer in polymer:
        natu_formatted_polymer.append(get_full_name(monomer))


    return natu_formatted_polymer

def convert_to_antismash_format(polymer: list[str]) -> list[str]:
    antismash_formatted_polymer = []
    for monomer in polymer:
        antismash_formatted_polymer.append(get_abbreviation(monomer))

    return antismash_formatted_polymer


def get_abbreviation(monomer: str) -> str:
    is_d = False
    has_nme = False
    has_me = False
    if monomer.startswith('D-') and 'D-lysergic acid' not in monomer:
        name = monomer[2:]
        is_d = True
    else:
        name = monomer

    if name.startswith('NMe-'):
        has_nme = True
        name = name[4:]
    if name.startswith('Me-'):
        has_me = True
        name = name[3:]

    try:
        short_name = get_substrate_by_name(name).short
    except ValueError:
        print(monomer)
        short_name = name

    if has_me:
        short_name = short_name # NATU currently does not handle C-methylations
    if has_nme:
        short_name = f"NMe-{short_name}"
    if is_d:
        short_name = f"D-{short_name}"
    return short_name


def get_full_name(monomer: str) -> str:
    is_d = False
    has_nme = False
    has_me = False
    if monomer.startswith('D-') and 'D-Lya' not in monomer:
        name = monomer[2:]
        is_d = True
    else:
        name = monomer

    if name.startswith('NMe-'):
        has_nme = True
        name = name[4:]
    if name.startswith('Me-'):
        has_me = True
        name = name[3:]

    try:
        long_name = get_substrate_by_name(name).long
    except ValueError:
        print(monomer)
        long_name = name

    if has_me:
        long_name = long_name # NATU currently does not handle C-methylations
    if has_nme:
        long_name = f"NMe-{long_name}"
    if is_d:
        long_name = f"D-{long_name}"
    return long_name
