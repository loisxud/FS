
beep <- function(n = 3){
  for(i in seq(n)){
    system("rundll32 user32.dll,MessageBeep -1")
    Sys.sleep(.5)
  }
}