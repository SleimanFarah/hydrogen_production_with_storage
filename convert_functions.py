def P_to_H2(P, LHV):  # P must be in kW
    H2 = P / LHV
    return H2


def H2_to_P(H2, LHV):  # H2 must be in kg
    P = LHV * H2
    return P