### find_activated_kinase.R ### 
# find activated kinase based on high phosphorylation level of downstream targets


### read in the kinase/substrate table/ phosphorylation data ### 

k_s_table = read.table()

pho_data = read.table()

# initiate an output table
out = c()

### loop through each of the kinase

unique_kinase = unique(t$kinase)

#find_substrate_pho = function(m1= , m2= )
for (kinase in unique_kinase){
  # find its substrate
  k_k_s_table = k_s_table[k_s_table$kinase == kinase,]
  substrates = k_k_s_table$substrate
  
  # get the phosphorylation data for the substrates
  substrates_pho_data = pho_data[row.names(pho_data) %in% substrates,]
  
  # get the summary statistic 
  
  for (sample in colnames(pho_data)){
    pho_data_s_s = pho_data[row.names(pho_data) %in% substratess,sample]
    # do comparison with the overall summary statistics
    
    # concatenate the result for that sample, kinase to the output table
    
  }
}
write.table(out, fn)