def fastcc(freq, alpha, detector=False, debug=False,dev=False):
  # Apply colour corrections to Planck LFI and WMAP data
  # freq should be one of 30, 44, 70, K, Ka, Q, V or W, as a string.
  # Actual frequencies are 28.4, 44.1 or 70.4 for LFI; 22.80 33.00 40.60 60.80 93.50 for WMAP
  # detector (optional) should be between '18' and '28' inclusive for LFI, as a string; they should of the style 'K11', 'K12', 'K1' for WMAP.
  # set 'debug' to see debug messages
  # set 'stable' to use the published WMAP and Planck numbers.
  # 
  # Version history:
  # Mike Peel   01-Feb-2013   v1.0 Initial version
  # Mike Peel   04-Feb-2013   v1.1 Correct convention
  # Locke Spencer 05-Feb-2013:  v1.2 Changed name to LFI_fastCC from planckcc (it only works for LFI)
  #                             Changed nested IF..ELSE groups to case statements
  #                             Changed invalid detector/frequency output to zero and removed internal stop in the code.
  # Mike Peel   24-Jul-2014   v2.0 Added WMAP (based on values from Paddy Leahy), renamed to fastcc. Note that function calls from v1 won't work with v2 without modification.
  # Mike Peel   13-Aug-2014   v2.1 Updated WMAP colour corrections, taking into account frequency drift during the mission (values from Paddy Leahy).
  # Mike Peel   06-Nov-2014   v2.2 Added 'dev' parameter (updated LFI and WMAP numbers from Paddy Leahy, including LFI bandpass shifts);
  #                                without this parameter the published values for LFI and WMAP will be returned.
  # Mike Peel   12-Mar-2015   v2.3 Update 70GHz values to fits to corrected colour corrections from Paddy Leahy
  # Mike Peel   23-Mar-2015   v2.4 Update WMAP colour corrections to corrected values from Paddy Leahy.
  # Mike Peel   22-Jan-2016   v2.5 Convert from IDL to Python

  # Define dictionaries containing the coefficients for different detectors and frequencies
  frequencies_v1 = {
    '30': [0.98520, 0.0131778, -0.00302],
    '44': [0.99059, 0.0079600, -0.00169],
    '70': [0.98149, 0.0152737, -0.00325],
    'K' : {'nu': [20.6, 22.8, 24.9], 'w': [0.332906, 0.374325, 0.292768], 'dT': 1.013438},
    'Ka': {'nu': [30.4, 33.0, 35.6], 'w': [0.322425, 0.387532, 0.290043], 'dT': 1.028413},
    'Q' : {'nu': [37.8, 40.7, 43.8], 'w': [0.353635, 0.342752, 0.303613], 'dT': 1.043500},
    'V' : {'nu': [55.7, 60.7, 66.2], 'w': [0.337805, 0.370797, 0.291399], 'dT': 1.098986},
    'W' : {'nu': [87.0, 93.5, 100.8], 'w': [0.337633, 0.367513, 0.294854], 'dT': 1.247521}
  }
  frequencies_v2 = {
    '30': [1.00513, 0.00301399, -0.00300699],
    '44': [0.994769, 0.00596703, -0.00173626],
    '70': [0.989711, 0.0106943, -0.00328671],
    'K' : [0.972902, 0.0190469, -0.00276464],
    'Ka': [0.983787, 0.0117567, -0.00183716],
    'Q' : [0.996854, 0.00496893, -0.00181359],
    'V' : [0.980322, 0.0143631, -0.00223596],
    'W' : [0.984848, 0.0112743, -0.00164595]
  }
  detectors_v1 = {
    '18': [0.98836, 0.0123556, -0.00394],
    '19': [0.93933, 0.0375844, -0.00225],
    '20': [0.95663, 0.0285644, -0.00273],
    '21': [0.97140, 0.0209690, -0.00318],
    '22': [1.02220,-0.0077263, -0.00327],
    '23': [1.00098, 0.0029940, -0.00240],
    '24': [0.99571, 0.0053247, -0.00175],
    '25': [0.98988, 0.0082248, -0.00161],
    '26': [0.98557, 0.0107023, -0.00175],
    '27': [0.98513, 0.0129780, -0.00288],
    '28': [0.98516, 0.0134605, -0.00318]
  }
  detectors_v2 = {
    '18': [0.977484, 0.0185055, -0.00391209],
    '19': [0.965314, 0.0234026, -0.00256943],
    '20': [0.968436, 0.0220869, -0.00285115],
    '21': [0.982854, 0.0142877, -0.00317682],
    '22': [1.049, -0.0237173, -0.00288312],
    '23': [0.990172, 0.0091968, -0.00238961],
    '1823': [0.983195, 0.0141778, -0.00317682],
    '1922': [1.00978, -0.000698302, -0.00328272],
    '2021': [0.97712, 0.0175904, -0.00308092],
    '24': [0.999958, 0.00309391, -0.00177223],
    '25': [0.994381, 0.00591109, -0.00162038],
    '26': [0.990046, 0.00854446, -0.00177223],
    '2526': [0.992115, 0.00717982, -0.00167233],
    '27': [1.00503, 0.00276424, -0.00282717],
    '28': [1.00491, 0.00334266, -0.00313287],
    'K11': [0.939366, 0.0346715, -0.00214346],
    'K12': [1.00894, 0.000982418, -0.00276923],
    'K1': [0.972902, 0.0190469, -0.00276464],
    'Ka11': [0.974784, 0.0159578, -0.00161958],
    'Ka12': [0.992978, 0.00737502, -0.00200839],
    'Ka1': [0.983787, 0.0117567, -0.00183716],
    'Q11': [0.990948, 0.00846474, -0.00204555],
    'Q12': [0.998159, 0.00404356, -0.00167233],
    'Q1': [0.994548, 0.00627672, -0.00186693],
    'Q21': [0.981607, 0.0126181, -0.00166893],
    'Q22': [1.01705, -0.00573297, -0.0016989],
    'Q2': [0.998986, 0.00378172, -0.00176723],
    'V11': [0.939474, 0.0354285, -0.00155105],
    'V12': [0.994737, 0.006396, -0.00217822],
    'V1': [0.966309, 0.0217416, -0.00209331],
    'V21': [1.00662, -0.000113686, -0.00217942],
    'V22': [0.977227, 0.0160255, -0.00220999],
    'V2': [0.991701, 0.0082012, -0.00226214],
    'W11': [0.988343, 0.00956424, -0.00211948],
    'W12': [0.9838, 0.0120015, -0.00173207],
    'W1': [0.986087, 0.0107974, -0.00193167],
    'W21': [0.978714, 0.0149705, -0.00148032],
    'W22': [0.992004, 0.00655744, -0.00146334],
    'W2': [0.985324, 0.0108262, -0.00149331],
    'W31': [0.977457, 0.0155997, -0.00131688],
    'W32': [0.993636, 0.0054001, -0.00134126],
    'W3': [0.985473, 0.0105855, -0.00135485],
    'W41': [0.991452, 0.0072962, -0.00181239],
    'W42': [0.973071, 0.0185705, -0.00153746],
    'W4': [0.982185, 0.0130277, -0.00170889]
  }

  # Pull out the desired coefficients from the above arrays
  if (detector != False):
    if (debug == True):
      print 'Using detector ',detector

    if (dev == True):
      cc = detectors_v2.get(detector, 0)
    else:
      cc = detectors_v1.get(detector, 0)

    if (cc == 0):
      print 'Invalid detector specified for fastcc, returning zero.'
      return 0

  else:
    if (debug == True):
      print 'Using frequency ',freq

    if (dev == True):
      cc = frequencies_v2.get(freq, 0)
    else:
      cc = frequencies_v1.get(freq, 0)

    if (cc == 0):
      print 'Invalid frequency specified for fastcc, returning zero.'
      return 0

  if (type(cc) is dict):
    # We have WMAP values
    beta=-2.0+alpha
    T0 = 1.0 * (cc['nu'][0]/cc['nu'][1]) ** beta
    T1 = 1.0
    T2 = 1.0 * (cc['nu'][2]/cc['nu'][1]) ** beta
    dT = 1.0 # Because conversion from T_CMB to T_RJ is done elsewhere.
    fastCC = 1.0 / (cc['dT'] * (cc['w'][0]*T0 + cc['w'][1]*T1 + cc['w'][2]*T2))
  else:
    fastCC = cc[0] + cc[1]*alpha + cc[2]*(alpha**2)

  return fastCC
