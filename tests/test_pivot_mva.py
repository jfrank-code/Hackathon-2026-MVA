import csv, subprocess, sys, tempfile, unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / 'src' / 'pivot_mva.py'

class PivotMVATest(unittest.TestCase):
    def test_bub1b_rank_and_complex_hla_suppression(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            pairs=root/'pairs.tsv'; hpo=root/'hpo.tsv'; inh=root/'inheritance.tsv'; out=root/'out'
            pairs.write_text(
                'PAIR_ID\tGENE\tCHROM_1\tPOS_1\tREF_1\tALT_1\tGT_1\tCONSEQUENCE_1\tCLINVAR_STATE_1\tCLINVAR_EXACT_1\tCLINVAR_DISEASES_1\tGNOMAD_AF_1\tDP_1\tGQ_1\tAF_READS_1\tCHROM_2\tPOS_2\tREF_2\tALT_2\tGT_2\tCONSEQUENCE_2\tCLINVAR_STATE_2\tCLINVAR_EXACT_2\tCLINVAR_DISEASES_2\tGNOMAD_AF_2\tDP_2\tGQ_2\tAF_READS_2\tINHERITANCE_MODEL_HINT\tINHERITANCE_HINT_REVIEWED\tINHERITANCE_CONDITION_HINT\tCLINVAR_PLP_COUNT\tCLINVAR_BENIGN_COUNT\tLOF_VARIANT_COUNT\tHPO_EXACT_MATCHES\tHPO_SPECIFIC_SUBTYPE_MATCHES\tHPO_ANY_RELATION_COUNT\tHPO_BEST_DISEASE_ID\n'
                'bub1b\tBUB1B\tchr15\t40209701\tT\tG\t0/1\tstop_gained\tPATHOGENIC_OR_LIKELY_PATHOGENIC\tYES\tMosaic_variegated_aneuploidy_syndrome_1\t7.8e-05\t46\t99\t0.54\tchr15\t40220612\tT\tG\t0/1\tmissense_variant\tNO_EXACT_ALLELE_MATCH_OR_UNCONFIRMED\tNO_OR_UNKNOWN\t\t6.8e-07\t28\t99\t0.46\tAR\tYES\tMosaic variegated aneuploidy syndrome 1\t1\t0\t1\t2\t1\t3\tOMIM:257300\n'
                'hla\tHLA-DQA1\tchr6\t32641443\tTGA\tT\t0/1\tframeshift_variant\tNO_EXACT_ALLELE_MATCH_OR_UNCONFIRMED\tNO_OR_UNKNOWN\t\tNA\t35\t99\t0.43\tchr6\t32641499\tTGG\tT\t0/1\tframeshift_variant\tNO_EXACT_ALLELE_MATCH_OR_UNCONFIRMED\tNO_OR_UNKNOWN\t\tNA\t37\t99\t0.51\tUNKNOWN\tNO\t\t0\t0\t2\t2\t0\t2\tOMIM:212750\n', encoding='utf-8')
            hpo.write_text('GENE\tDISEASE_ID\tN_EXACT\tN_SPECIFIC_SUBTYPE_RELATIONS\tN_TERMS_WITH_ANY_ONTOLOGY_RELATION\tMATCHED_PHENOTYPES\nBUB1B\tOMIM:257300\t2\t1\t3\tRhabdomyosarcoma=MORE_SPECIFIC_SUBTYPE; Short stature=EXACT; Small for gestational age=EXACT\nHLA-DQA1\tOMIM:212750\t2\t0\t2\tShort stature=EXACT; Failure to thrive=EXACT\n', encoding='utf-8')
            inh.write_text('GENE\tCONDITION\tINHERITANCE\tSOURCE_OR_NOTE\tREVIEWED\nBUB1B\tMVA1\tAR\tClinGen\tYES\n', encoding='utf-8')
            p=subprocess.run([sys.executable,str(SCRIPT),'--pairs',str(pairs),'--hpo-summary',str(hpo),'--inheritance',str(inh),'--out',str(out)],capture_output=True,text=True)
            self.assertEqual(p.returncode,0,msg=p.stdout+'\n'+p.stderr)
            with (out/'pivot_mva_scored_pairs.tsv').open(encoding='utf-8',newline='') as f:
                scored=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(scored[0]['GENE'],'BUB1B')
            self.assertEqual(scored[0]['PIVOT_ASSOCIATION_CLASS'],'MENDELIAN_DISEASE')
            self.assertEqual(scored[0]['PIVOT_INHERITANCE_STATUS'],'REVIEWED_AR_TWO_HET_PHASE_UNKNOWN')
            hla=next(x for x in scored if x['GENE']=='HLA-DQA1')
            self.assertEqual(hla['PIVOT_HPO_DOMAIN_COUNT'],'0')
            self.assertEqual(hla['PIVOT_AF_KNOWN'],'0') # missing AF is not zero-valued frequency
            with (out/'pivot_mva_submission_primary_only.csv').open(encoding='utf-8',newline='') as f:
                draft=list(csv.DictReader(f))
            self.assertEqual(len(draft),1)
            self.assertEqual(draft[0]['gene'] if 'gene' in draft[0] else draft[0]['pos_1'],'40209701')
            self.assertGreater(float(draft[0]['epcr']),0)

if __name__=='__main__': unittest.main()
