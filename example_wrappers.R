metric <- function(x){
  1/x
}
metric(0)

metric_wrapper <- function(x){
  if(x == 0){
    x <- x + 1e-06
  }
  metric(x)
}
metric_wrapper(0)

metric_wrapper <- function(x){
  tryCatch({
    if(x == 0) stop()
    metric(x)
  }, error = function(e){ NA })
}
metric_wrapper(0)
