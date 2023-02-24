



# TODO
def json2df(json_in):
    """
    Parse json file. (It really depends on how deep the json file is).
    Args:
        json_in: str
            input json format file.

    Returns:
        return a string but formatted like a dataframe.
    """
    import json
    fh = open(json_in)
    # with open(json_in) as fp:
    # out.write('Head line.\n')
    json_file = json.load(fh)
    # out_df = pd.DataFrame(index=None, columns=None)
    out_st_df = ''
    n0 = json_file['name']
    for cld1 in json_file['children']:
        n1 = cld1['name']
        for cld2 in cld1['children']:
            n2 = cld2['name']
            for cld3 in cld2['children']:
                n3 = cld3['name']  # "ko01000" is strange. what this for?
                if 'children' not in cld3:
                    # out.write(f'{n0}\t{n1}\t{n2}\t{n3}\t{n4}\n')
                    # out.write(f'{n0}\t{n1}\t{n2}\t{n3}\n')
                    # out_df.append(pd.DataFrame([n0, n1, n2, n3]).T)
                    out_st_df += f'{n0}\t{n1}\t{n2}\t{n3}\n'
                else:
                    # print(cld3)
                    for cld4 in cld3['children']:
                        n4 = cld4['name']
                        if 'children' not in cld4:
                            out_st_df += f'{n0}\t{n1}\t{n2}\t{n3}\t{n4}\n'
                        else:
                            for cld5 in cld4['children']:
                                n5 = cld5['name']
                                if 'children' not in cld5:
                                    out_st_df += f'{n0}\t{n1}\t{n2}\t{n3}\t{n4}\t{n5}\n'
                                else:
                                    for cld6 in cld5['children']:
                                        n6 = cld6['name']
                                        out_st_df += f'{n0}\t{n1}\t{n2}\t{n3}\t{n4}\t{n5}\t{n6}\n'
    return out_st_df
