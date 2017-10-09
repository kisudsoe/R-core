---
Title: RSQLite.md
Date: 2017-10-08 SUN
Version: v1.0
Note : Here is example code for SQLite connection in R
---

```r
#install.packages("RSQLite",repos='http://cran.us.r-project.org')
install.packages("sqldf",repos='http://cran.us.r-project.org')

# Sample db examples
library(sqldf)
db = datasetsDb(); print(dbListTables(db))
head(dbReadTable(db,"CO2")); print(dim(dbReadTable(db,"CO2")))
head(dbGetQuery(db,"SELECT * FROM CO2")); print(dim(dbGetQuery(db,"SELECT * FROM CO2")))

# Make sqlite db from csv file
read.csv.sql(file="iris.csv",
             sql=c("attach 'iris.sqlite' as new", # create *.sqlite file
                   "CREATE TABLE irisdb AS SELECT * FROM file"),
             dbname="iris.sqlite")
print(dbListTables(db)) # Check db list
sqldf("select * from sqldb limit5", dbname="iris.sqlite")

# Make sqlite db from data.frame variable
iris2 = read.csv("iris.csv"); print(dim(iris2))
db = dbConnect(SQLite(), dbname="iris.sqlite")
#dbSendQuery(conn=db,"create table irisdb (sepal_l,sepal_w,petal_l,petal_w)")
dbSendQuery(conn=db,"drop table if exists irisdb") # Erase previous table
dbWriteTable(db, "irisdb", iris2)
print(dbListTables(db)) # Check db list
print(dbListFields(db,'irisdb')) # Column in table
sqldf("SELECT * FROM irisdb limit 5", dbname="iris.sqlite") # Table head
#dbReadTable(db,'irisdb2') # Table values

dbDisconnect(db) # Close connection
```
