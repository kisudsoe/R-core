# 190808 - 1st generation for Mito-GTEx project
# This function is for collect DAVID GO result *.txt files
# and save as *.csv or *.rds file.
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
readgo = function(
	dir = "",
	select_filt = "", #"*_CC.txt"
	filt_term = "", #
	write_csv = T,
	save_rds = T
) {
	# Get file paths
    paths = list.files(dir); #length(paths)
	paths_select = grep(select_filt,paths,value=T); #print(paths_select)
    
    # Read DAVID GO txt files
	txt_li = lapply(paths_select,function(p){
		f_p = paste0(dir,'/',p)
		name = unlist(strsplit(p,'_')) %>%
			head(-2) %>% paste0(collapse='_')
			#%>% print
		dat = read.delim(f_p)
		dat$Organs = rep(name,nrow(dat))
		return(dat)
	})
	txt_df = ldply(txt_li,data.frame) # plyr package
	#txt_df %>% dim %>% print
    
    # Filt and save GO contents
	if(filt_term!="") txt_df = filter(txt_df,grepl("mitochond",Term))
	#dim(txt_df) %>% print; unique(txt_df$Term) %>% print
	out = txt_df$Organs %>% table %>% data.frame
	if(write_csv==T) {
		f_name1 = paste0(dir,'_',filt_term,'.csv')
		write.csv(txt_df,f_name1)
		cat(paste0('>> file generate: ',f_name1,'\n'))
	} 
	if(save_rds==T) {
		f_name2 = paste0(dir,'_',filt_term,'.rds')
		saveRDS(txt_df,f_name2)
		cat(paste0('>> file generate: ',f_name2,'\n'))
	}
    return(out)
}