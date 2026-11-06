"""Explicit deterministic output extension: lumped reaction and frozen relaxation.

Generated H derivatives and physical history updates remain unchanged. Previous
nodal damage is the generator's U-DU field; gamma uses accepted GP H only.
No abaqus_ufl package/template changes are made.
"""
def transform(src):
    raw_lines = len(src.splitlines())
    replaced_lines = replacement_lines = added_lines = 0

    def replace(old, new, *, addition=False):
        nonlocal src, replaced_lines, replacement_lines, added_lines
        assert src.count(old) == 1, (old, src.count(old))
        src = src.replace(old, new)
        if addition:
            added_lines += len(new.splitlines())-len(old.splitlines())
        else:
            replaced_lines += len(old.splitlines())
            replacement_lines += len(new.splitlines())

    assert 'DOUBLE PRECISION :: d_old(4)' in src
    assert 'd_old(ii_v) = U(idx) - DU(idx, 1)' in src
    replace('      DOUBLE PRECISION :: phase_storage_real',
            '      DOUBLE PRECISION :: phase_storage_real\n'
            '      DOUBLE PRECISION :: phase_relaxation', addition=True)
    replace('        phase_storage_real = DBLE(phase_storagez_eval)',
            '        phase_storage_real = DBLE(phase_storagez_eval)\n'
            '        phase_relaxation = PROPS(13)*MIN(hhist_old,PROPS(14))',
            addition=True)
    replace('     &      - phase_storage_real*sh4(ii_v)*wdetJ',
            '     &      - (d_node(ii_v)-(1.0d0-PROPS(5))*vfrac_node(ii_v)\n'
            '     &      * phase_storage_real+phase_relaxation\n'
            '     &      * (d_node(ii_v)-d_old(ii_v)))*sh4(ii_v)*wdetJ')
    for suffix in ['_dF(k,l)', '_dc', '_dd', '_dvfrac']:
        replace('+ dphase_storage'+suffix+'\n',
                '- (1.0d0-PROPS(5))*vfrac_node(ii_v)\n'
                '     &        * dphase_storage'+suffix+'\n')
    anchor = 'C       drp_resid_dF (sign=+1, value-grad)'
    diagonal = (
        'C       Exact nodal derivatives; relaxation is frozen accepted-history data\n'
        '        DO ii_v = 1, 4\n'
        '          row = edof_d(ii_v)\n'
        '          col = edof_d(ii_v)\n'
        '          AMATRX(row,col) = AMATRX(row,col)\n'
        '     &      +(1.0d0+phase_relaxation)*sh4(ii_v)*wdetJ\n'
        '          col = edof_vfrac(ii_v)\n'
        '          AMATRX(row,col) = AMATRX(row,col)\n'
        '     &      -(1.0d0-PROPS(5))*phase_storage_real*sh4(ii_v)*wdetJ\n'
        '        END DO\n\n')
    replace(anchor, diagonal+anchor, addition=True)
    assert len(src.splitlines()) == raw_lines-replaced_lines+replacement_lines+added_lines
    return src, dict(
        description='Lump nodal d/vfrac reaction with generated GP H chain rule; '
                    'add frozen gamma*(d_i-dold_i), gamma=eta*min(Hold,Hcap), '
                    'and exact damage diagonal',
        original_lines_replaced=replaced_lines, replacement_lines=replacement_lines,
        additional_lines=added_lines)
