import awkward as ak

def get_trigger_mask(events, is_MC=False, is_UL=False):
            # MC: always pass trigger
            if is_MC:
                return ak.ones_like(events.event, dtype=bool)
        
            # UL data
            if is_UL:
                return (
                    events.HLT.AK8PFJet400_TrimMass30 |
                    events.HLT.AK8PFJet420_TrimMass30 |
                   events.HLT.PFHT1050 |
                    events.HLT.PFJet500
                )
        
            # non-UL (Run 2 / pre-UL)
            return (
                events.HLT.AK8PFJetFwd40 |
                events.HLT.AK8PFJet230_SoftDropMass40 |
                events.HLT.QuadPFJet103_88_75_15 |
                events.HLT.QuadPFJet105_88_76_15 |
                events.HLT.QuadPFJet111_90_80_15 |
                events.HLT.DiPFJetAve260 |
                events.HLT.PFJetFwd400 |
                events.HLT.AK8PFJet420_MassSD30 |
                events.HLT.AK8DiPFJet270_270_MassSD30 |
                events.HLT.AK8PFJet450_MassSD30 |
                events.HLT.AK8PFJet450_SoftDropMass40 |
                events.HLT.AK8PFJet425_SoftDropMass40 |
                events.HLT.AK8PFJet500
            )
