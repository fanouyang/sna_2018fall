class_data<- read.csv("/Users/fanouyang/Desktop/CI5301_2015Fall_analysis/class_data/orig_class_dis.csv")
class_data_new<-class_data %>% select(vert1_id,vert2_id)
class_data_new<- read.csv("/Users/fanouyang/Desktop/CI5301_2015Fall_analysis/class_data/class_data_new.csv")
head(class_data_new)
#write.csv(class_data_new,file = "/Users/fanouyang/Desktop/CI5301_2015Fall_analysis/class_data/class_data_new.csv")

class_data_new %>% group_by(vert1_id,vert2_id) %>% summarize(fre=n())

library(reshape2)
class_data_new<-class_data_new %>% dcast(vert1_id~vert2_id)

#write.csv(class_data_new,file = "/Users/fanouyang/Desktop/CI5301_2015Fall_analysis/class_data/class_data_new.csv")



class(class_data_new)
class_net=network(class_data_new)
class_net


gplot(class_net,gmode="graph")
