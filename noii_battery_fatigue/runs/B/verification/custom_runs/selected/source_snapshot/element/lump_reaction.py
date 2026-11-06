"""Deterministic generator-output extension for nodal reaction quadrature.

Only damage reaction assembly changes. H and its derivatives remain generated;
H/alpha/W history updates, all other equations and gradient quadrature are intact.
No abaqus_ufl package/template modification is made.
"""
def transform(src):
    original = "     &      - phase_storage_real*sh4(ii_v)*wdetJ"
    replacement = ("     &      - (d_node(ii_v)-(1.0d0-PROPS(5))*vfrac_node(ii_v)\n"
                   "     &      * phase_storage_real)*sh4(ii_v)*wdetJ")
    assert src.count(original) == 1
    src = src.replace(original, replacement)
    for suffix in ['_dF(k,l)', '_dc', '_dd', '_dvfrac']:
        token = '+ dphase_storage'+suffix+'\n'
        assert src.count(token) == 1, (suffix, src.count(token))
        src = src.replace(token, '- (1.0d0-PROPS(5))*vfrac_node(ii_v)\n'
                                '     &        * dphase_storage'+suffix+'\n')
    anchor = 'C       drp_resid_dF (sign=+1, value-grad)'
    assert src.count(anchor) == 1
    diagonal = ("C       Exact nodal partial derivatives of lumped phase reaction\n"
                "        DO ii_v = 1, 4\n"
                "          row = edof_d(ii_v)\n"
                "          col = edof_d(ii_v)\n"
                "          AMATRX(row,col) = AMATRX(row,col)+sh4(ii_v)*wdetJ\n"
                "          col = edof_vfrac(ii_v)\n"
                "          AMATRX(row,col) = AMATRX(row,col)\n"
                "     &      -(1.0d0-PROPS(5))*phase_storage_real*sh4(ii_v)*wdetJ\n"
                "        END DO\n\n")
    src = src.replace(anchor, diagonal+anchor)
    return src, ('Lump phase reaction to nodal d/vfrac with GP H; multiply generated '
                 'H chain-rule tangent by -a*vfrac_node(row) and add nodal diagonals')
