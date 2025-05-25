\ Model matlab
\ LP format - for model browsing. Use MPS format to capture full model detail.
\ Signature: 0x7f4d3cb40fe8155e
Minimize
  10 x1 + 7 x2 + 16 x3 + 6 x4 + 40 x5 + 24 x6 + 4 x7 + 45 x8 + 27 x9
   + 4.5 x10 + 32 x11 + 19.2 x12 + 3.2 x13 + 55 x14 + 33 x15 + 5.5 x16
   + 1000 x17 + 1000 x18 + 1000 x19 + 1000 x20
Subject To
 R0: x1 + x2 + x3 + x4 >= 15
 R1: 10 x1 + 7 x2 + 16 x3 + 6 x4 <= 220
 R2: - x1 + x5 + x6 + x7 - x17 <= 0
 R3: - x2 + x8 + x9 + x10 - x18 <= 0
 R4: - x3 + x11 + x12 + x13 - x19 <= 0
 R5: - x4 + x14 + x15 + x16 - x20 <= 0
 R6: x5 + x8 + x11 + x14 >= 5
 R7: x6 + x9 + x12 + x15 >= 4.000025
 R8: x7 + x10 + x13 + x16 >= 3.001325
Bounds
End
