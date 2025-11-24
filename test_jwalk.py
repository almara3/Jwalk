from Jwalk.jwalk import runJwalk


def test_runJwalk():
    max_dist = 60
    vox = 1
    surface = True
    xl_path = 'Jwalk/Examples/1FGA_xl_list.txt'
    pdb_path = 'Jwalk/Examples/1FGA.pdb'

    results = runJwalk(pdb_path, max_dist, vox, surface, xl_path, ncpus=1)

    print(results)


if __name__ == "__main__":
    test_runJwalk()
