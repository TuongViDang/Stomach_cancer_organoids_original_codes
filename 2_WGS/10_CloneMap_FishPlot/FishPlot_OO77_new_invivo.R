#FishPlots don't need to change clone name

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(fishplot)
source("/group/poetsch_projects/poetsch_sc/Conifer/script/fishplot_package/object.R")

CONIPHER_DIR="/group/poetsch_projects/poetsch_sc/Conifer/tree30_change_name"
plot_dir=file.path(CONIPHER_DIR,"FishPlot_new_invivo")
dir.create(plot_dir, recursive = T)

patient = "OO77"

tree_df = read.delim(file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient,"output/tree_30_final/Trees/allTrees.txt"),comment.char = "#",header=F)
colnames(tree_df) = c("parent","clusterID")
tree_df = rbind(c(0,10),tree_df)

Find_direct_children <- function(cluster,tree){
    direct_children = tree %>% filter(parent == cluster) %>% pull(clusterID)
    direct_children
}

Find_all_children <- function(cluster,tree) {
  all_children <- vector()
  direct_children <- Find_direct_children(cluster,tree=tree)
  all_children <- c(all_children, direct_children)
  for (child in direct_children) {
    grandchildren <- Find_all_children(child, tree = tree)
    all_children <- c(all_children, grandchildren)
  }
  return(all_children)
}

Clones_all = read.delim(file.path("/group/poetsch_projects/poetsch_sc/Conifer",patient,"output/tree_30_final/Trees/cloneProportionsMinErrorTrees.txt"))
Clones_all <- Clones_all  %>% left_join(tree_df,by = "clusterID")


New_proportion <- function(cluster,sample,clones_df,tree){
    children = Find_all_children(cluster,tree = tree)
    children = c(cluster,children)
    children_in_sample = clones_df[(clones_df[sample] > 0 ) & clones_df$clusterID %in% children,]
    Sum_Pro_in_Sample = sum(children_in_sample[sample])
    return(Sum_Pro_in_Sample)
    
}

tumor_new = purrr::map(Clones_all$clusterID,~ New_proportion(.x,sample = paste0(patient,"_tumor"),clones_df = Clones_all, tree=tree_df)) %>% unlist()
ex_new = purrr::map(Clones_all$clusterID,~ New_proportion(.x,sample = paste0(patient,"_ex"),clones_df = Clones_all, tree=tree_df)) %>% unlist()
in_new = purrr::map(Clones_all$clusterID,~ New_proportion(.x,sample = paste0(patient,"_in"),clones_df = Clones_all, tree=tree_df)) %>% unlist()

Clones_all_new = data.frame(clusterID = Clones_all$clusterID, parent= Clones_all$parent,tumor = tumor_new, ex_vivo = ex_new, in_vivo = in_new)

#rename the cluster
clusterID <- Clones_all_new$clusterID
names(clusterID) <- rownames(Clones_all_new)
parentIDs_new = c()
 for (x in Clones_all_new$parent){ 
    if (x == 0){x_new = 0}
    else{
     x_new = names(clusterID)[clusterID == x] }
    parentIDs_new = c(parentIDs_new,x_new)
  }
cols <-  c("darkkhaki", "dodgerblue4", "lightgreen", "mediumseagreen", "lightpink",  "lightgoldenrod1","cadetblue2","plum") 

names(cols) <- clusterID

#Untreated 
    timepoints=c(0,1)      


    frac.table = as.matrix(Clones_all_new[,c(3,3)])
    #provide a vector listing each clone's parent
    #(0 indicates no parent)
    parents = as.integer(parentIDs_new)

    #create a fish object
    fish = createFishObject(frac.table,parents,timepoints=timepoints)

    #calculate the layout of the drawing
    fish = layoutClones(fish)
    

    #draw the plot, using the splining method (recommended)
    #and providing both timepoints to label and a plot title

    fish@col = cols


    pdf(file.path(plot_dir,paste("fish_plot_",patient,"_Begin_Untreated.pdf")),width= 8, height = 6.5)
    fishPlot(fish,shape="spline",title.btm=patient, pad.left = 4, ramp.angle = 0.7,
             cex.title=0.7, 
             vlab=c( "",""),bg.col = c("white","white","white"))

    dev.off()

#Untreated and ex vivo
    timepoints=c(0,1)      


    frac.table = as.matrix(Clones_all_new[,c(3,4)])
    #provide a vector listing each clone's parent
    #(0 indicates no parent)
    parents = as.integer(parentIDs_new)

    #create a fish object
    fish = createFishObject(frac.table,parents,timepoints=timepoints)

    #calculate the layout of the drawing
    fish = layoutClones(fish)
    

    #draw the plot, using the splining method (recommended)
    #and providing both timepoints to label and a plot title

    names(cols) <- clusterID
    fish@col = cols


    pdf(file.path(plot_dir,paste("fish_plot_",patient,"_EX.pdf")),width= 13.5, height = 6.5)
    fishPlot(fish,shape="spline",title.btm=patient, pad.left = 1, ramp.angle = 3,
             cex.title=0.7, vlines=c(0,1), 
            # vlab=c( "Untreated","Ex vivo treated"),
            bg.col = c("white","white","white"))

    dev.off()
    
    

#In vivo (new as Anna's suggestion)
    timepoints=c(0,1)      
    
    #create new column as pseudo in vivo only contain clone 2 (clone 3 in this table)
    Clones_all_new$pseudo_in = c(0,0,100,0,0,0,0,0)
    frac.table = as.matrix(Clones_all_new[,c(6,6)])

    #provide a vector listing each clone's parent(0 indicates no parent)
    parents = as.integer(parentIDs_new)
    #change the parent of 3 to 0
    parents[3] = 0
    #create a fish object
    fish = createFishObject(frac.table,parents,timepoints=timepoints)

    #calculate the layout of the drawing
    fish = layoutClones(fish)
    fish@col = cols


    pdf(file.path(plot_dir,paste0("fish_plot_",patient,"_IN.pdf")),width= 13.5, height = 6.5)
    fishPlot(fish,shape="spline",title.btm=patient, pad.left = 1, ramp.angle = 3,
             cex.title= 0.7,
            #vlines=c(0,1), 
            # vlab=c( "Untreated","In vivo treated"),
            bg.col = c("white","white","white"))

    dev.off()
