def dimension_cusp_forms_sp4z(k, J):
    """
    Uses Tsushima's formula from Theorem 4 of 'An explicit dimension formula
    for the spaces of generalized automorphic forms with respect to Sp(2,Z)'
    Proc. Japan Acad. 59 Ser A (1983).

    Tsushima proves the correctness of his formula for (j = 0 and  k >= 4)
    or (j >= 1 and k >= 5), but Bergstroem-Faber-van der Geer prove that it
    holds for (j >= 0 and k >= 4), see page 97 of 'Siegel modular forms of
    degree three and the cohomology of local systems' Sel. Math. New Ser.
    (2014) 20:83-124.
    """
    if (J % 2) == 1:
        return ZZ(0)
    j = ZZ(J/2)
    if j < 0:
        raise ValueError("j cannot be negative")

    k = ZZ(k)
    if k < 4:
        raise ValueError("not implemented for k < 4")

    res = 2^(-7) * 3^(-3) * 5^(-1) * (2*j+1) * (k-2) * (2*j+k-1) * (2*j+2*k-3)
    res += - 2^(-5) * 3^(-2) * (2*j+1) * (2*j+2*k-3)
    res += 2^(-4) * 3^(-1) * (2*j+1)
    res += (-1)^k * (2^(-7) * 3^(-2) * 7 * (k-2) * (2*j+k-1) - 2^(-4) * 3^(-1) * (2*j+2*k-3) + 2^(-5) * 3)
    res += (-1)^j * (2^(-7) * 3^(-1) * 5 * (2*j+2*k-3) - 2^(-3))
    res += (-1)^k * (-1)^j * 2^(-7) * (2*j+1)

    i = CyclotomicField(4).gen()
    rho = CyclotomicField(3).gen()
    omega = CyclotomicField(5).gen()
    sigma = CyclotomicField(12).gen()

    res += (i^k * (2^(-6) * 3^(-1) * i * (2*j+k-1) - 2^(-4) * i)).trace()
    res += ((-1)^k * i^j * 2^(-5) * (i+1)).trace()
    res += (i^k * (-1)^j * (2^(-6) * 3^(-1) * (k-2) - 2^(-4))).trace()
    res += ((-i)^k * i^j * 2^(-5) * (i+1)).trace()
    res += ((-1)^k * rho^j * 3^(-3) * (rho+1)).trace()
    res += (rho^k * rho^j * 2^(-2) * 3^(-4) * (2*rho+1) * (2*j+1)).trace()
    res += - (rho^k * (-rho)^j * 2^(-2) * 3^(-2) * (2*rho+1)).trace()
    res += ((-rho)^k * rho^j * 3^(-3)).trace()
    res += (rho^j * (2^(-1) * 3^(-4) * (1-rho) * (2*j+2*k-3) - 2^(-1) * 3^(-2) * (1-rho))).trace()
    res += (rho^k * (2^(-3) * 3^(-4) * (rho+2) * (2*j+k-1) - 2^(-2) * 3^(-3) * (5*rho+6))).trace()
    res += - ((-rho)^k * (2^(-3) * 3^(-3) * (rho+2) * (2*j+k-1) - 2^(-2) * 3^(-2) * (rho+2))).trace()
    res += (rho^k * (rho^2)^j * (2^(-3) * 3^(-4) * (1-rho) * (k-2) + 2^(-2) * 3^(-3) * (rho-5))).trace()
    res += ((-rho)^k * (rho^2)^j * (2^(-3) * 3^(-3) * (1-rho) * (k-2) - 2^(-2) * 3^(-2) * (1-rho))).trace()
    res += (omega^k * (omega^4)^j * 5^(-2)).trace()
    res += - (omega^k * (omega^3)^j * 5^(-2) * omega^2).trace()
    res += ((sigma^7)^k * (-1)^j * 2^(-3) * 3^(-2) * (sigma^2+1)).trace()
    res += - ((sigma^7)^k * (sigma^8)^j * 2^(-3) * 3^(-2) * (sigma+sigma^3)).trace()
    return ZZ(res)


def generating_series_cusp_forms_sp4z_wt3(j):
    """
    From page 45 of Petersen 'Cohomology of
    local systems on the moduli of principally polarized abelian surfaces'
    """
    R.<x> = PowerSeriesRing(QQ, default_prec=j+1)
    num = x^36
    denom = (1-x^6) * (1-x^8) * (1-x^10) * (1-x^12)
    f = num / denom
    return f


def generating_series_modular_forms_sp4z_wt4(j):
    """
    From page 8 of Ibukiyama 'Lifting conjectures from vector valued Siegel
    modular forms of degree two'
    """
    R.<x> = PowerSeriesRing(QQ, default_prec=j+1)
    num = x^8+x^12-x^18-x^20-x^22+x^28+x^30+x^32
    denom = (1-x^6) * (1-x^8) * (1-x^10) * (1-x^12)
    f = num / denom
    return f


def generating_series_cusp_forms_sp4z_wt4(j):
    """
    From page 8 of Ibukiyama 'Lifting conjectures from vector valued Siegel
    modular forms of degree two'
    """
    R.<x> = PowerSeriesRing(QQ, default_prec=j+1)
    num = x^24 * (1+x^4+x^8-x^10)
    denom = (1-x^6) * (1-x^8) * (1-x^10) * (1-x^12)
    f = num / denom
    return f


def generating_series_cusp_forms_sp4z_wt5(j):
    """
    From page 115 of Ibukiyama "A conjecture on a Shimura type correspondence
    for Siegel modular forms, and Harder's conjecture on congruences"
    """
    R.<x> = PowerSeriesRing(QQ, default_prec=j+1)
    num = x^18+x^20+x^24
    denom = (1-x^6) * (1-x^8) * (1-x^10) * (1-x^12)
    f = num / denom
    return f


def generating_series_cusp_forms_sp4z_wt7(j):
    """
    From page 115 of Ibukiyama "A conjecture on a Shimura type correspondence
    for Siegel modular forms, and Harder's conjecture on congruences"
    """
    R.<x> = PowerSeriesRing(QQ, default_prec=j+1)
    num = x^12+x^14+x^16+x^18+x^20
    denom = (1-x^6) * (1-x^8) * (1-x^10) * (1-x^12)
    f = num / denom
    return f


def generating_series_cusp_forms_sp4z_odd_wt(k, j):
    """
    From page 115 of Ibukiyama "A conjecture on a Shimura type correspondence
    for Siegel modular forms, and Harder's conjecture on congruences"

    Note that there is a shift in the parameters j and k in Ibukiyama's formula
    """
    R.<x,y> = PowerSeriesRing(QQ, default_prec=k+j+2)
    denom = (1-y^2) * (1-y^3) * (1-y^5) * (1-y^6) * (1-x^3) * (1-x^4) * (1-x^5) * (1-x^6)
    # this needs work, the numerator seems to only be given for a particular
    # range in Ibukiyama
    f = num / denom
    return f


def generating_series_modular_forms_sp4z_j10_even_wt(k):
    """
    From Lemma 7.1 of Takemori 'Structure theorems for vector valued Siegel
    modular forms of degree 2 and weight det^k otimes Sym(10)'
    """
    R.<y> = PowerSeriesRing(QQ, default_prec=k+1)
    num = y^6+y^8+2*y^10+2*y^12+3*y^14+2*y^16+y^18+y^20-y^24-y^26
    denom = (1-y^4) * (1-y^6) * (1-y^10) * (1-y^12)
    f = num / denom
    return f


def generating_series_modular_forms_sp4z_j10_odd_wt(k):
    """
    From Lemma 7.1 of Takemori 'Structure theorems for vector valued Siegel
    modular forms of degree 2 and weight det^k otimes Sym(10)'
    """
    R.<y> = PowerSeriesRing(QQ, default_prec=k+1)
    num = y^9+y^11+y^13+3*y^15+3*y^17+2*y^19+y^21+y^23-y^27-y^29
    denom = (1-y^4) * (1-y^6) * (1-y^10) * (1-y^12)
    f = num / denom
    return f


def generating_series_modular_forms_sp4z_j0_even_wt(k):
    """
    From Section 9 of van der Geer 'Siegel modular forms'
    """
    R.<y> = PowerSeriesRing(QQ, default_prec=k+1)
    num = 1
    denom = (1-y^4) * (1-y^6) * (1-y^10) * (1-y^12)
    f = num / denom
    return f


def generating_series_modular_forms_sp4z_j0_odd_wt(k):
    """
    From Section 9 of van der Geer 'Siegel modular forms'
    """
    R.<y> = PowerSeriesRing(QQ, default_prec=k+1)
    num = y^35
    denom = (1-y^4) * (1-y^6) * (1-y^10) * (1-y^12)
    f = num / denom
    return f


def generating_function_numerator_cusp_forms_sp4z_k(k, prec=None):
    denom = generating_function_denominator_cusp_forms_sp4z_k(k)
    if prec is None:
        prec = denom.degree()
    lst = [dimension_cusp_forms_sp4z(k, j) for j in range(prec)]
    R.<x> = PowerSeriesRing(QQ)
    res = (R(lst) * denom).add_bigoh(prec)
    return res.polynomial().change_ring(ZZ)


def generating_function_denominator_cusp_forms_sp4z_k(k):
    R.<x> = PolynomialRing(ZZ)
    res = (1-x^6) * (1-x^8) * (1-x^10) * (1-x^12)
    return res


def generating_function_numerator_cusp_forms_sp4z_j(j, prec=None):
    denom = generating_function_denominator_cusp_forms_sp4z_j(j)
    if prec is None:
        prec = denom.degree()
    lst = [dimension_cusp_forms_sp4z(k, j) for k in range(prec)]
    R.<x> = PowerSeriesRing(QQ)
    res = (R(lst) * denom).add_bigoh(prec)
    return res.polynomial().change_ring(ZZ)


def generating_function_denominator_cusp_forms_sp4z_j(j):
    R.<x> = PolynomialRing(ZZ)
    res = (1-x^4) * (1-x^6) * (1-x^10) * (1-x^12)
    return res


def dimension_modular_forms_sp4z(k, j):
    """
    """
    j = ZZ(j)

    if j < 0:
        raise ValueError("j cannot be negative")

    k = ZZ(k)
    if k < 0 or (k == 0 and j > 0):
        return ZZ(0)

    if (k % 2) == 1 or (k == 2):
        return dimension_cusp_forms_sp4z(k, j)

    # if we're here then k is even and k >= 4
    res = dimension_cusp_forms_sp4z(k, j) + dimension_cusp_forms(1, k+j)
    if j == 0:
        # there's one Siegel Eisenstein form
        res += 1
    return res


def mods(ts, ell, m):
    return ts[(m % ell)]


def dimension_cusp_forms_Gamma_e(k, j):
    """
    From Theorem 6.2 in Ibukiyama and Wakatsuki
    """
    H1e = 2^(-6) * 3^(-3) * 5^(-1) * (j+1) * (k-2) * (j+k-1) * (j+2*k-3)
    H1u = - 2^(-6) * 3^(-2) * (j+1) * (j+2*k-3) + 2^(-4) * 3^(-1) * (j+1)
    H1 = H1e + H1u

    H2e = 2^(-6) * 3^(-2) * (-1)^k * (j+k-1) * (k-2)
    H2qu = - 2^(-4) * 3^(-1) * (-1)^k * (j+2*k-3) + 2^(-6) * 3 * (-1)^k
    H2 = H2e + H2qu

    H3e = 0
    H3qu = - 2^(-3) * mods([(-1)^(j/2), -1, (-1)^(j/2+1), 1], 4, k) + 2^(-4) * mods([1, (-1)^(j/2), -1, (-1)^(j/2+1)], 4, k)
    H3 = H3e + H3qu

    H4e = 2^(-2) * 3^(-3) * (mods([j+k-1, -(j+k-1), 0], 3, k) + mods([k-2, 0, -(k-2)], 3, j+k))
    H4qu = - 2^(-2) * 3^(-2) * (mods([1, -1, 0], 3, k) + mods([1, 0, -1], 3, j+k)) - 3^(-2) * (mods([0, -1, -1], 3, k) + mods([1, 1, 0], 3, j+k))
    H4 = H4e + H4qu

    H5e = 2^(-2) * 3^(-2) * (mods([-(j+k-1), -(j+k-1), 0, j+k-1, j+k-1, 0], 6, k) + mods([k-2, 0, -(k-2), -(k-2), 0, k-2], 6, j+k))
    H5qu = - 2^(-2) * 3^(-1) * (mods([-1, -1, 0, 1, 1, 0], 6, k) + mods([1, 0, -1, -1, 0, 1], 6, j+k))
    H5 = H5e + H5qu

    H6e = 2^(-6) * (-1)^(j/2) * (j+2*k-3) + 2^(-6) * (-1)^(j/2+k) * (j+1)
    H6qu = - 2^(-3) * (-1)^(j/2)
    H6 = H6e + H6qu

    H7e = 3^(-3) * (j+2*k-3) * mods([1, -1, 0], 3, j) + 2^(-1) * 3^(-3) * (j+1) * mods([0, 1, -1], 3, j+2*k)
    H7qu = - 2^(-1) * 3^(-1) * mods([1, -1, 0], 3, j)
    H7 = H7e + H7qu

    H8 = 0

    if (j % 6) == 0:
        h9 = mods([1, 0, 0, -1, 0, 0], 6, k)
    elif (j % 6) == 2:
        h9 = mods([-1, 1, 0, 1, -1, 0], 6, k)
    elif (j % 6) == 4:
        h9 = mods([0, -1, 0, 0, 1, 0], 6, k)
    H9 = 2^(-1) * 3^(-2) * h9

    if (j % 10) == 0:
        h10 = mods([1, 0, 0, -1, 0], 5, k)
    elif (j % 10) == 2:
        h10 = mods([-1, 1, 0, 0, 0], 5, k)
    elif (j % 10) == 4:
        h10 = 0
    elif (j % 10) == 6:
        h10 = mods([0, 0, 0, 1, -1], 5, k)
    elif (j % 10) == 8:
        h10 = mods([0, -1, 0, 0, 1], 5, k)
    H10 = 2 * 5^(-1) * h10

    if (j % 8) == 0:
        h11 = mods([1, 1, -1, -1], 4, k)
    elif (j % 8) == 2:
        h11 = mods([-1, 1, 1, -1], 4, k)
    elif (j % 8) == 4:
        h11 = mods([-1, -1, 1, 1], 4, k)
    elif (j % 8) == 6:
        h11 = mods([1, -1, -1, 1], 4, k)
    H11 = 2^(-3) * h11

    H12 = 0

    res = H1 + H2 + H3 + H4 + H5 + H6 + H7 + H8 + H9 + H10 + H11 + H12

    return res


def dimension_cusp_forms_sp4z_sgn(k, j):
    d1 = dimension_cusp_forms_Gamma_e(k, j)
    d2 = dimension_cusp_forms_sp4z(k, j)
    return d1 - d2


def dimension_V(k, j):
    if (k % 2) == 0:
        d1 = 2^(-1) * dimension_cusp_forms(1, k+j/2)
        d2 = 2^(-1) * sum([dimension_cusp_forms(1, k+j-2*a) * dimension_cusp_forms(1, k+2*a) for a in range(j/2+1)])

    if (k % 2) == 1:
        d1 = - 2^(-1) * dimension_cusp_forms(1, k+j/2)
        d2 = 2^(-1) * sum([dimension_cusp_forms(1, k-1+j-2*a) * dimension_cusp_forms(1, k+1+2*a) for a in range(j/2)])

    return d1 + d2


def dimension_V2(k, j):
    d1 = 2^(-1) * (-1)^k * dimension_cusp_forms(1, k+j/2)
    d2 = 2^(-1) * sum([dimension_cusp_forms(1, k+j-a) * dimension_cusp_forms(1, k+a) for a in range(j+1)])
    return d1 + d2


def dimension_W(k, j):
    d1 = 2^(-1) * (-1)^(k+1) * dimension_modular_forms(1, k+j/2-6)
    d2 = 2^(-1) * sum([dimension_modular_forms(1, k+j-a-6) * dimension_modular_forms(1, k+a-6) for a in range(j+1)])
    return d1 + d2
