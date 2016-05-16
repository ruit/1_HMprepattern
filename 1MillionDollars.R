# Tian R. <ruitian@yeah.net>
# May 16, 2016
# 1 Million dollars behind 3 doors problem


#Stick to the initial choice the winning odds is 1/3, whereras by switching you get 2/3

#based on the R simulation, it is seen that by switching you get twice the winning odds.


threeDoor<-function(N){

	#label the 3 doors
	doors<-c("a", "b", "c")
	
	result<-c()

	for (i in 1:N){
		#trial times N

		money<-sample(doors)[1] #sample the door, "a" or "b" or "c" behind which there is the 1 million dollars

		guess<-sample(doors)[1] #sample the first choice
	
		#"open"	the door which neither chosen nor with 1 million dollars behind, this is critical!	
		open<-sample(doors[which(doors != guess & doors != money)])[1]

		#switch to the door which is not the chosen one or the opened one
		switch.opt<-doors[which(doors  != guess & doors !=open)]


		#summarize the winning results
		if(guess==money){result<-c(result, "stayWin")} #no switching
	
		if(switch.opt==money){result<-c(result,"switchWin")} #switching
	}

	return (list(stayWinRate=length(which(result=="stayWin"))/N, switchWinRate=length(which(result=="switchWin"))/N))
}


#simulate 1000 times

threeDoor(1000)

# a result
#$stayWinRate
#[1] 0.339

#$swithchWinRate
#[1] 0.661


