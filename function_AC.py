

def Annualized_cost(LoanRate, LifeExpectancy, Cost, Maintenance):
    if LoanRate !=0:
        CRF=LoanRate/(1-(1+LoanRate)**(-LifeExpectancy))
    else:
        CRF=1

    AC=(CRF*Cost+Maintenance*LifeExpectancy)/LifeExpectancy
    return AC
