
#Use the libraries
library(jsonlite) #https://cran.r-project.org/web/packages/jsonlite/
library(httr)

#Fill in the AphiaID you need
AphiaID <- 127160

#Build the URL to get the data from
url <- sprintf("https://www.marinespecies.org/rest/AphiaClassificationByAphiaID/%d", AphiaID);

#Get the actual data from the URL
classificationTree <- fromJSON(url)

#Walk the classification tree
currentTreeItem = classificationTree
while (!is.null(currentTreeItem )) {
  print(sprintf("ID: %d, RANK: %s, ScientificName: %s",
                currentTreeItem$AphiaID,
                currentTreeItem$rank,
                currentTreeItem$scientificname
  ));
  #You can access the variables individually
  #print(currentTreeItem$AphiaID);
  #print(currentTreeItem$scientificname);
  #print(currentTreeItem$rank);

  #Get next item in the tree
  currentTreeItem <- currentTreeItem$child;
}
