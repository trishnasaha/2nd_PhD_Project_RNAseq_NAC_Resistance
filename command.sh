## file processing for sorting out drug target for hub gene

~/Dropbox/Trishna_Projects/03_Text_Mining/analysis_of_genes/01_Updated_Oct_29/Drug_Repurposing$ while read id; do echo $id ; grep -w $id repurposing_drugs_20200324.txt ; done < hub_genes.txt > repurposing_drugs_20200324_hub_genes.txt  from the repurposing drug dataset


## file processing for sorting out drug target for hub gene from the drug dataset

~/Dropbox/Trishna_Projects/03_Text_Mining/analysis_of_genes/Drug/TTD_db$ grep -e "TARGETID" -e "GENENAME" -e "DRUGINFO" P1-01-TTD_target_download.txt | grep -v -e "Phase" -e "Investigative" -e "Terminated" -e "Patented" -e "Withdrawn" -e "#NAME?" -e "Clinical trial" -e "Application submitted" -e "Registered" -e "Discontinued" | grep -v -e "^TARGETID" -e "^GENENAME" -e "^DRUGINFO" | cut -f2,3,4,5 | sed 's/TARGETID/\nTARGETID/g' > P1-01-TTD_target_download.tsv 

trishna@trishna-linux:~/Dropbox/Trishna_Projects/04_Drug_designing/02_drug_vina$ grep "^@<TRIPOS>MOLECULE" fda.mol2 | wc -l
2115

obabel -imol2 fda.mol2 -opdbqt -O fda.pdbqt

trishna@trishna-linux:~/Dropbox/Trishna_Projects/04_Drug_designing/02_drug_vina$ grep -c "^MODEL" fda.pdbqt 
2115

trishna@trishna-linux:~/Dropbox/Trishna_Projects/04_Drug_designing/02_drug_vina$ vina_split --input fda.pdbqt 
Prefix for ligands will be fda_ligand_
Prefix for flexible side chains will be fda_flex_
