### validate_kinase_sub.R ### 
# look at correlations of kinase and downstream substrates phosphorylation status
# sub-P ~ a*kinase + k
# sub-P ~ a*kinase + b*sub + k

### read in the kinase/substrate table/ phosphorylation data ### 
K_S_f = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/Phospho_databases/PhosphositePlus/data/Kinase_Substrate_Dataset_human_final_hugoified.txt"
k_s_table = read.table()
BRCA_pro_f = ""
BRCA_pho_f = "/Users/khuang/Box\ Sync/PhD/proteogenomics/CPTAC_pan3Cancer/pan3can_shared_data/BRCA/TCGA_Breast_BI_Phosphoproteome.phosphosite.itraq_abbrev_normlized.tsv"

pho_data = read.table("")

# initiate an output table
out = c()

### loop through each of the kinase

unique_kinase = unique(t$kinase)

results;
#find_substrate_pho = function(m1= , m2= )
for (kinase in unique_kinase){
  # find its substrate
  k_k_s_table = k_s_table[k_s_table$kinase == kinase,]
  substrates = k_k_s_table$substrate
  
  for (substrate in substrates){
    # fit a regression model 
    fit = glm()
    results = fdfjd
    fit2 = fjdl
    results2 = fdjla
    
    if ( P < 0.05){
      plot()
    }
  }
  
}

write.table(out, fn)

# plotting module