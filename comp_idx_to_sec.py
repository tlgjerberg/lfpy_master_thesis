
def idx_to_sec_conversion(num_key='0'):
    """
    Convert comprtment index to segment name
    """
    idx_to_sec = {
        '0': 'soma',
        '471': 'apic',
        '104': 'axon',
        '125': 'axon_term',
        '42': 'axon_term2'
    }

    return idx_to_sec.get(num_key, "no index")
