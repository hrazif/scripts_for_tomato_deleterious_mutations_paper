EA00585_x_het=read.table("EA00585_x_het.txt", header=T, check.names=F, stringsAsFactors = F)


EA00585_x_het$ANC="none"
EA00585_x_het$DER="none"

temp1=EA00585_x_het[EA00585_x_het$EA00585=="0/0",]
temp2=EA00585_x_het[EA00585_x_het$EA00585=="1/1",]


EA00585_x_het[EA00585_x_het$EA00585=="0/0","ANC"]=temp1$REF
EA00585_x_het[EA00585_x_het$EA00585=="0/0","DER"]=temp1$ALT

EA00585_x_het[EA00585_x_het$EA00585=="1/1","ANC"]=temp2$ALT
EA00585_x_het[EA00585_x_het$EA00585=="1/1","DER"]=temp2$REF



write.table(EA00585_x_het, "EA00585_x_het_wANC_DER.txt", sep="\t", quote=F, row.names = F, col.names = T)