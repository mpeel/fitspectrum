# Test code for fastcc
# 
# Version history:
# Mike Peel   01-Feb-2013   v1.0 Initial version
# Mike Peel   04-Feb-2013   v1.1 Update format
# Locke Spencer 05-Feb-2013:  v1.2 changed planckcc to LFI_fastcc within this code to follow changes to LFI_fastcc
#                             renamed this routine LFI_fastcc_test, from planckcc_test, to follow convention of other changes.
#                             changed this routine from a script to a procedure so that it can be included within the hfi_lfi_test_script example routine also.
# Mike Peel   24-Jul-2014   v2.0 Expand to include WMAP, and to use new function calls.
# Mike Peel   06-Nov-2014   v2.1 Add development version.
# Mike Peel   22-Jan-2016   v2.2 Convert from IDL to Python. Does not currently work on arrays of spectra.

from fastcc import fastcc

spectra = -2.0
#spectra = [-2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

print 'STABLE VERSION:'
print 'Detector	alpha',spectra
print 'LFI-18',fastcc('70',spectra,detector='18')
print 'LFI-19',fastcc('70',spectra,detector='19')
print 'LFI-20',fastcc('70',spectra,detector='20')
print 'LFI-21',fastcc('70',spectra,detector='21')
print 'LFI-22',fastcc('70',spectra,detector='22')
print 'LFI-23',fastcc('70',spectra,detector='23')
print '70GHz',fastcc('70',spectra)

print 'LFI-24',fastcc('44',spectra,detector='24')
print 'LFI-25',fastcc('44',spectra,detector='25')
print 'LFI-26',fastcc('44',spectra,detector='26')
print '44GHz',fastcc('44',spectra)

print 'LFI-27',fastcc('30',spectra,detector='27')
print 'LFI-28',fastcc('30',spectra,detector='28')
print '30GHz',fastcc('30',spectra)

print 'WMAP'
print 'K',fastcc('K',spectra)
print 'Ka',fastcc('Ka',spectra)
print 'Q',fastcc('Q',spectra)
print 'V',fastcc('V',spectra)
print 'W',fastcc('W',spectra)

print 'DEVELOPMENT VERSION:'
print 'Detector	alpha',spectra
print 'LFI-18',fastcc('70',spectra,detector='18',dev=True)
print 'LFI-19',fastcc('70',spectra,detector='19',dev=True)
print 'LFI-20',fastcc('70',spectra,detector='20',dev=True)
print 'LFI-21',fastcc('70',spectra,detector='21',dev=True)
print 'LFI-22',fastcc('70',spectra,detector='22',dev=True)
print 'LFI-23',fastcc('70',spectra,detector='23',dev=True)
print 'LFI-18-23',fastcc('70',spectra,detector='1823',dev=True)
print 'LFI-19-22',fastcc('70',spectra,detector='1922',dev=True)
print 'LFI-20-21',fastcc('70',spectra,detector='2021',dev=True)
print '70GHz',fastcc('70',spectra,dev=True)

print 'LFI-24',fastcc('44',spectra,detector='24',dev=True)
print 'LFI-25',fastcc('44',spectra,detector='25',dev=True)
print 'LFI-26',fastcc('44',spectra,detector='26',dev=True)
print 'LFI-25-26',fastcc('44',spectra,detector='2526',dev=True)
print '44GHz',fastcc('44',spectra,dev=True)

print 'LFI-27',fastcc('30',spectra,detector='27',dev=True)
print 'LFI-28',fastcc('30',spectra,detector='28',dev=True)
print '30GHz',fastcc('30',spectra,dev=True)

print 'WMAP'
print 'K11',fastcc('K',spectra,detector='K11',dev=True)
print 'K12',fastcc('K',spectra,detector='K12',dev=True)
print 'K1',fastcc('K',spectra,detector='K1',dev=True)
print 'K',fastcc('K',spectra,dev=True)

print 'Ka11',fastcc('Ka',spectra,detector='Ka11',dev=True)
print 'Ka12',fastcc('Ka',spectra,detector='Ka12',dev=True)
print 'Ka1',fastcc('Ka',spectra,detector='Ka1',dev=True)
print 'Ka',fastcc('Ka',spectra,dev=True)

print 'Q11',fastcc('Q',spectra,detector='Q11',dev=True)
print 'Q12',fastcc('Q',spectra,detector='Q12',dev=True)
print 'Q1',fastcc('Q',spectra,detector='Q1',dev=True)
print 'Q21',fastcc('Q',spectra,detector='Q21',dev=True)
print 'Q22',fastcc('Q',spectra,detector='Q22',dev=True)
print 'Q2',fastcc('Q',spectra,detector='Q2',dev=True)
print 'Q',fastcc('Q',spectra,dev=True)

print 'V11',fastcc('V',spectra,detector='V11',dev=True)
print 'V12',fastcc('V',spectra,detector='V12',dev=True)
print 'V1',fastcc('V',spectra,detector='V1',dev=True)
print 'V21',fastcc('V',spectra,detector='V21',dev=True)
print 'V22',fastcc('V',spectra,detector='V22',dev=True)
print 'V2',fastcc('V',spectra,detector='V2',dev=True)
print 'V',fastcc('V',spectra,dev=True)

print 'W11',fastcc('W',spectra,detector='W11',dev=True)
print 'W12',fastcc('W',spectra,detector='W12',dev=True)
print 'W1',fastcc('W',spectra,detector='W1',dev=True)
print 'W21',fastcc('W',spectra,detector='W21',dev=True)
print 'W22',fastcc('W',spectra,detector='W22',dev=True)
print 'W2',fastcc('W',spectra,detector='W2',dev=True)
print 'W31',fastcc('W',spectra,detector='W31',dev=True)
print 'W32',fastcc('W',spectra,detector='W32',dev=True)
print 'W3',fastcc('W',spectra,detector='W3',dev=True)
print 'W41',fastcc('W',spectra,detector='W41',dev=True)
print 'W42',fastcc('W',spectra,detector='W42',dev=True)
print 'W4',fastcc('W',spectra,detector='W4',dev=True)
print 'W',fastcc('W',spectra,dev=True)
