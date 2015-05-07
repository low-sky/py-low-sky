# This is the suite of function that checks whether
#  we're improving agreement between data and model
# and penalizes moves we don't like.  

def Objective(data, rms, model):
    chisq = (data-model)**2/(2*rms**2)
# Blanked data?
    return chisq.sum()
