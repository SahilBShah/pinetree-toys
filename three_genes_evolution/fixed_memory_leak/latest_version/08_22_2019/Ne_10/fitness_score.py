import math

#Calculates the fitness of the new mutation
def calc_fitness(variant_fit, orig_fit):
    """
    Calculates the fitness of the new mutation and compares it to the fitness of the old mutation.
    The 'orig_fit' is the old fitness value previously calculated.
    The 'variant_fit' is the new fitness value calculated to determine if mutation is accepted or rejected.
    """

    xi = orig_fit
    xj = variant_fit
    Ne = 10

    if xj == xi:
        return(1.0 / float(Ne))
    if xj > xi:
        return(1.0)
    else:
        variant_fit = (xj / xi) ** (2 * Ne - 2)
        return variant_fit

    try:
        resolved = ((1-pow((xi/xj), 2)) / (1-pow((xi/xj), (2 * float(Ne)))))
    except OverflowError as e:
        resolved = 0.0
    return (resolved)
