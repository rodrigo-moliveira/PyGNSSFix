# TODO: add the appropriate unit test setup...

from src.data_types.date import Epoch
from src.data_types.date.frames import get_gmst


def validation_test():
    print("Test to validate Greenwich Mean Sidereal Time (GMST) computation against official IAU SOFA code.")
    epoch_test = Epoch(2007, 4, 5, 15, 5, 23, scale="UTC")
    epoch_ut1 = epoch_test.change_scale("UT1")
    epoch_tt = epoch_test.change_scale("TT")
    gmst82 = get_gmst(epoch_test, 'IAU82')
    gmst00 = get_gmst(epoch_test, 'IAU00')
    gmst06 = get_gmst(epoch_test, 'IAU06')

    # IAU REFS
    IAU_JD_UTC = 2454195.500000000000000+0.628738425925926
    IAU_JD_UT1 = 2454195.500000000000000+0.628737599274306
    IAU_JD_TT = 2454195.500000000000000+0.629492870370370
    IAU_GMST82 = 1.045176704334423
    IAU_GMST00 = 1.045176678373470
    IAU_GMST06 = 1.045176677938862

    print(f"JD UTC: computed value = {epoch_test.jd} | IAU SOFA = {IAU_JD_UTC} -> diff = {epoch_test.jd - IAU_JD_UTC}")
    print(f"JD UT1: computed value = {epoch_ut1.jd} | IAU SOFA = {IAU_JD_UT1} -> diff = {epoch_ut1.jd - IAU_JD_UT1}")
    print(f"JD TT: computed value = {epoch_tt.jd} | IAU SOFA = {IAU_JD_TT} -> diff = {epoch_tt.jd - IAU_JD_TT}")
    print(f"gmst82: computed value = {gmst82} | IAU SOFA = {IAU_GMST82} -> diff = {gmst82 - IAU_GMST82}")
    print(f"gmst00: computed value = {gmst00} | IAU SOFA = {IAU_GMST00} -> diff = {gmst00 - IAU_GMST00}")
    print(f"gmst06: computed value = {gmst06} | IAU SOFA = {IAU_GMST06} -> diff = {gmst06 - IAU_GMST06}")


if __name__ == "__main__":
    validation_test()
