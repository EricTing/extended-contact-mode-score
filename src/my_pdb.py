from Bio.PDB.PDBIO import Select


def buildSelect(selected_residues):
    chain_id = selected_residues[0][0]
    residue_ids = [_[-1] for _ in selected_residues]

    class MySelect(Select):
        def accept_chain(self, chain):
            if chain.get_id() == chain_id:
                return True
            else:
                return False

        def accept_residue(self, residue):
            if residue.get_id()[1] in residue_ids:
                return True
            else:
                return False

    return MySelect
