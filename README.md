# Asian-Option-pricing
Asian option pricer for an option that pays the maximum average price out of three correlated stocks minus the strike price.
Price determined using Monte carlo simulation. The price paths to maturity are simulated for the 3 stocks, and for
each iteration of the monte carlo simulation a payoff for the option is calculated. Finally the average of the payoffs is
calculated and discounted to today. 
