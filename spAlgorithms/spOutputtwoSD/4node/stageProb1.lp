\ENCODING=ISO-8859-1
\Problem name: 4node_1

Minimize
 obj1: D0.A.B + D0.B.A + D1.A.B + D1.A.E + D1.B.E + D2.A.B + D2.A.C + D2.B.C
       + D3.A.C + D3.C.A + D4.A.C + D4.A.E + D4.C.E + D5.A.C + D5.A.B + D5.C.B
       + D6.B.A + D6.A.B + D7.B.A + D7.B.C + D7.A.C + D8.B.C + D8.B.A + D8.C.A
       + D9.B.C + D9.C.B + D10.B.C + D10.B.E + D10.C.E + D11.B.E + D11.E.B
       + D12.B.E + D12.B.C + D12.E.C + D13.E.C + D13.C.E + D14.E.C + D14.E.B
       + D14.C.B + D15.E.C + D15.E.A + D15.C.A + D16.E.B + D16.B.E + D17.E.B
       + D17.E.C + D17.B.C + D18.E.B + D18.E.A + D18.B.A + D19.C.A + D19.A.C
       + D20.C.A + D20.C.B + D20.A.B + D21.C.B + D21.B.C + D22.C.B + D22.C.A
       + D22.B.A + D23.C.B + D23.C.E + D23.B.E + D24.C.E + D24.E.C + D25.C.E
       + D25.C.B + D25.E.B + T0.ABC + T0.ABE + T1.ABC + T2.ABE + T3.ACB
       + T3.ACE + T4.ACB + T5.ACE + T6.BAC + T8.BCE + T9.BCA + T9.BCE + T10.BCA
       + T11.BEC + T13.ECA + T13.ECB + T14.ECA + T15.ECB + T16.EBA + T16.EBC
       + T17.EBA + T18.EBC + T19.CAB + T21.CBA + T21.CBE + T22.CBE + T23.CBA
       + T24.CEB + S0.B.A + S1.B.E + S2.B.C + S3.C.A + S4.C.E + S5.C.B + S6.A.B
       + S7.A.C + S8.C.A + S9.C.B + S10.C.E + S11.E.B + S12.E.C + S13.C.E
       + S14.C.B + S15.C.A + S16.B.E + S17.B.C + S18.B.A + S19.A.C + S20.A.B
       + S21.B.C + S22.B.A + S23.B.E + S24.E.C + S25.E.B + 1300 YAB + 1300 YAC
       + 1300 YAE + 1300 YBA + 1300 YBC + 1300 YBE + 1300 YCA + 1300 YCB
       + 1300 YCE + 1300 YEA + 1300 YEB + 1300 YEC
Subject To
 CBAB:  D0.A.B + D1.A.B + D2.A.B + D5.A.B + D6.A.B + D20.A.B + T3.ACB + T4.ACB
        + YAB >= 5
 CBAC:  D2.A.C + D3.A.C + D4.A.C + D5.A.C + D7.A.C + D19.A.C + T0.ABC + T1.ABC
        + YAC >= 8
 CBAE:  D1.A.E + D4.A.E + T0.ABE + T2.ABE + T3.ACE + T5.ACE + YAE >= 4
 CBBA:  D0.B.A + D6.B.A + D7.B.A + D8.B.A + D18.B.A + D22.B.A + T9.BCA
        + T10.BCA + YBA >= 4.5
 CBBC:  D2.B.C + D7.B.C + D8.B.C + D9.B.C + D10.B.C + D12.B.C + D17.B.C
        + D21.B.C + T6.BAC + T11.BEC + YBC >= 6.7
 CBBE:  D1.B.E + D10.B.E + D11.B.E + D12.B.E + D16.B.E + D23.B.E + T8.BCE
        + T9.BCE + YBE >= 4.2
 CBCA:  D3.C.A + D8.C.A + D15.C.A + D19.C.A + D20.C.A + D22.C.A + T21.CBA
        + T23.CBA + YCA >= 10.1
 CBCB:  D5.C.B + D9.C.B + D14.C.B + D20.C.B + D21.C.B + D22.C.B + D23.C.B
        + D25.C.B + T19.CAB + T24.CEB + YCB >= 8.3
 CBCE:  D4.C.E + D10.C.E + D13.C.E + D23.C.E + D24.C.E + D25.C.E + T21.CBE
        + T22.CBE + YCE >= 12.2
 CBEA:  D15.E.A + D18.E.A + T13.ECA + T14.ECA + T16.EBA + T17.EBA + YEA >= 3.2
 CBEB:  D11.E.B + D14.E.B + D16.E.B + D17.E.B + D18.E.B + D25.E.B + T13.ECB
        + T15.ECB + YEB >= 5.3
 CBEC:  D12.E.C + D13.E.C + D14.E.C + D15.E.C + D17.E.C + D24.E.C + T16.EBC
        + T18.EBC + YEC >= 7.5
 TBAB:  T19.CAB - S6.A.B - S20.A.B  = 0
 TBAC:  T6.BAC - S7.A.C - S19.A.C  = 0
 TBBA:  T16.EBA + T17.EBA + T21.CBA + T23.CBA - S0.B.A - S18.B.A - S22.B.A  = 0
 TBBC:  T0.ABC + T1.ABC + T16.EBC + T18.EBC - S2.B.C - S17.B.C - S21.B.C  = 0
 TBBE:  T0.ABE + T2.ABE + T21.CBE + T22.CBE - S1.B.E - S16.B.E - S23.B.E  = 0
 TBCA:  T9.BCA + T10.BCA + T13.ECA + T14.ECA - S3.C.A - S8.C.A - S15.C.A  = 0
 TBCB:  T3.ACB + T4.ACB + T13.ECB + T15.ECB - S5.C.B - S9.C.B - S14.C.B  = 0
 TBCE:  T3.ACE + T5.ACE + T8.BCE + T9.BCE - S4.C.E - S10.C.E - S13.C.E  = 0
 TBEB:  T24.CEB - S11.E.B - S25.E.B  = 0
 TBEC:  T11.BEC - S12.E.C - S24.E.C  = 0
 PB0A:  D0.A.B + T0.ABC + T0.ABE + Z0A  = 0
 PB1A:  D1.A.B + D1.A.E + T1.ABC + Z1A  = 0
 PB2A:  D2.A.B + D2.A.C + T2.ABE + Z2A  = 0
 PB3A:  D3.A.C + T3.ACB + T3.ACE + Z3A  = 0
 PB4A:  D4.A.C + D4.A.E + T4.ACB + Z4A  = 0
 PB5A:  D5.A.C + D5.A.B + T5.ACE + Z5A  = 0
 PB6B:  D6.B.A + T6.BAC + Z6B  = 0
 PB7B:  D7.B.A + D7.B.C + Z7B  = 0
 PB8B:  D8.B.C + D8.B.A + T8.BCE + Z8B  = 0
 PB9B:  D9.B.C + T9.BCA + T9.BCE + Z9B  = 0
 PB10B: D10.B.C + D10.B.E + T10.BCA + Z10B  = 0
 PB11B: D11.B.E + T11.BEC + Z11B  = 0
 PB12B: D12.B.E + D12.B.C + Z12B  = 0
 PB13E: D13.E.C + T13.ECA + T13.ECB + Z13E  = 0
 PB14E: D14.E.C + D14.E.B + T14.ECA + Z14E  = 0
 PB15E: D15.E.C + D15.E.A + T15.ECB + Z15E  = 0
 PB16E: D16.E.B + T16.EBA + T16.EBC + Z16E  = 0
 PB17E: D17.E.B + D17.E.C + T17.EBA + Z17E  = 0
 PB18E: D18.E.B + D18.E.A + T18.EBC + Z18E  = 0
 PB19C: D19.C.A + T19.CAB + Z19C  = 0
 PB20C: D20.C.A + D20.C.B + Z20C  = 0
 PB21C: D21.C.B + T21.CBA + T21.CBE + Z21C  = 0
 PB22C: D22.C.B + D22.C.A + T22.CBE + Z22C  = 0
 PB23C: D23.C.B + D23.C.E + T23.CBA + Z23C  = 0
 PB24C: D24.C.E + T24.CEB + Z24C  = 0
 PB25C: D25.C.E + D25.C.B + Z25C  = 0
 PB0B:  - D0.A.B + D0.B.A - T0.ABC - T0.ABE + S0.B.A - Z0A + Z0B  = 0
 PB1B:  - D1.A.B + D1.B.E - T1.ABC + S1.B.E - Z1A + Z1B  = 0
 PB2B:  - D2.A.B + D2.B.C - T2.ABE + S2.B.C - Z2A + Z2B  = 0
 PB3C:  - D3.A.C + D3.C.A - T3.ACB - T3.ACE + S3.C.A - Z3A + Z3C  = 0
 PB4C:  - D4.A.C + D4.C.E - T4.ACB + S4.C.E - Z4A + Z4C  = 0
 PB5C:  - D5.A.C + D5.C.B - T5.ACE + S5.C.B - Z5A + Z5C  = 0
 PB6A:  - D6.B.A + D6.A.B - T6.BAC + S6.A.B - Z6B + Z6A  = 0
 PB7A:  - D7.B.A + D7.A.C + S7.A.C - Z7B + Z7A  = 0
 PB8C:  - D8.B.C + D8.C.A - T8.BCE + S8.C.A - Z8B + Z8C  = 0
 PB9C:  - D9.B.C + D9.C.B - T9.BCA - T9.BCE + S9.C.B - Z9B + Z9C  = 0
 PB10C: - D10.B.C + D10.C.E - T10.BCA + S10.C.E - Z10B + Z10C  = 0
 PB11E: - D11.B.E + D11.E.B - T11.BEC + S11.E.B - Z11B + Z11E  = 0
 PB12E: - D12.B.E + D12.E.C + S12.E.C - Z12B + Z12E  = 0
 PB13C: - D13.E.C + D13.C.E - T13.ECA - T13.ECB + S13.C.E - Z13E + Z13C  = 0
 PB14C: - D14.E.C + D14.C.B - T14.ECA + S14.C.B - Z14E + Z14C  = 0
 PB15C: - D15.E.C + D15.C.A - T15.ECB + S15.C.A - Z15E + Z15C  = 0
 PB16B: - D16.E.B + D16.B.E - T16.EBA - T16.EBC + S16.B.E - Z16E + Z16B  = 0
 PB17B: - D17.E.B + D17.B.C - T17.EBA + S17.B.C - Z17E + Z17B  = 0
 PB18B: - D18.E.B + D18.B.A - T18.EBC + S18.B.A - Z18E + Z18B  = 0
 PB19A: - D19.C.A + D19.A.C - T19.CAB + S19.A.C - Z19C + Z19A  = 0
 PB20A: - D20.C.A + D20.A.B + S20.A.B - Z20C + Z20A  = 0
 PB21B: - D21.C.B + D21.B.C - T21.CBA - T21.CBE + S21.B.C - Z21C + Z21B  = 0
 PB22B: - D22.C.B + D22.B.A - T22.CBE + S22.B.A - Z22C + Z22B  = 0
 PB23B: - D23.C.B + D23.B.E - T23.CBA + S23.B.E - Z23C + Z23B  = 0
 PB24E: - D24.C.E + D24.E.C - T24.CEB + S24.E.C - Z24C + Z24E  = 0
 PB25E: - D25.C.E + D25.E.B + S25.E.B - Z25C + Z25E  = 0
End
