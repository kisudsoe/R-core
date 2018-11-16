# pdtime function

`11/16/2018` - first write

This function is to calculate processed time and to put time stamp of running code. You need to add `t = Sys.time()` at top of your code, and call pdtime function at bottom.

```R
pdtime = function(time) {
    t=Sys.time()
    d=difftime(t,time,unit="sec")
    if(d<60) d=paste0(round(d,1)," sec, ",t)
    else if(d>=60 && d<3600) d=paste0(round(d/60,1)," min, ",t)
    else if(d>=3600) d=paste0(round(d/3600,1)," hr, ",t)
    return(d)
}
```

