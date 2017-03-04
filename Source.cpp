#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace std;


double maxi(double a, double b, double c);		//Returns the maximum out of 3 numbers
double * twoSNormals(int flag);					//Generates independent pairs of standard normal random variable
double cumulative_normal(double x);				//Approximation to the cumulative normal function (Hastings) used for the black Scholes formula
double black_scholes(int type_of_option, double S0, double K, double T, double r, double vol);	//Returns the Black Scholes price of a call or put option

class asianOptionPricing {

private:

	double S0_stock1;			//Initial Stock 1 price
	double S0_stock2;			//Initial Stock 2 price
	double S0_stock3;			//Initial Stock 3 price
	double K;					//Strike
	double T;					//Maturity
	double r;					//Risk free rate
	double vol1;				//Stock 1 volatility
	double vol2;				//Stock 2 volatility
	double vol3;				//Stock 3 volatility
	double cor12;				//Correlation between Stock 1 and 2
	double cor13;				//Correlation between Stock 1 and 3
	double cor23;				//Correlation between Stock 2 and 3
	double past_days_average;	//Number of days from maturity to calculate average

	void houseKeeping();


public:

	int type; //Type is used in the main function and thus needs to be declared in the public scope of the class

	//Constructor with parameters
	asianOptionPricing(int typ,double S0_1, double S0_2, double S0_3, double strike,
		double maturity, double risk_free, double vola1, double vola2, double vola3,
		double corel12, double corel13, double corel23, double av_days) {

		type = typ;
		S0_stock1 = S0_1;
		S0_stock2 = S0_2;
		S0_stock3 = S0_3;
		K = strike;
		T = maturity;
		r = risk_free;
		vol1 = vola1;
		vol2 = vola2;
		vol3 = vola3;
		cor12 = corel12;
		cor13 = corel13;
		cor23 = corel23;
		past_days_average = av_days;

		//Function to fix any wrong input
		houseKeeping();

	}

	//Constructor for user input
	asianOptionPricing() {
		
		cout << "Please enter the number 1 for call or 2 for put: " << endl;
		cin >> type;
		cout << "Enter the 3 current stock prices (separated by space): ";
		cin >> S0_stock1 >> S0_stock2 >> S0_stock3;
		cout << "Enter volatilities of each stock: ";
		cin >> vol1 >> vol2 >> vol3;
		cout << "Enter strike price, maturity and risk free rate: ";
		cin >> K >> T >> r;
		cout << "Enter correlations (1-2, 1-3, 2-3): ";
		cin >> cor12 >> cor13 >> cor23;
		cout << "Number of days for the average price from maturity: ";
		cin >> past_days_average;

		//Function to fix any wrong input
		houseKeeping();
	}

	void priceOption();

};

void asianOptionPricing::houseKeeping(){

	//Variable type must have value 1 or 2, otherwise ask user for re-input
	while (type < 1 || type > 2){ 
		cout << "Re-enter call or put: " << endl;
		cin >> type;

	}
	
	//Correlations should be between -1 and 1
	while (fabs(cor12) > 1 || fabs(cor13) > 1 || fabs(cor23) > 1) {

		cout << "Re-enter correlations (1-2, 1-3, 2-3): ";
		cin >> cor12 >> cor13 >> cor23;

	}

	while (past_days_average > T * 252 || past_days_average < 0)	{
		cout << "Re-enter Number of days for the average price from maturity : ";
		cin >> past_days_average;
	}

	//This can be done for all the parameters of the option but these are the most important
	//For instance the prices, strike, risk free and volatilities cant be less than 0
}

void asianOptionPricing::priceOption() {

	//initialize variables
	double opt_price = 0, opt_price_antithetic = 0;
	int const n = 10000; //Number of simulations for the Monte Carlo
	int i, j;
	double average1 = 0, average2 = 0, average3 = 0;
	double average1_anti = 0, average2_anti = 0, average3_anti = 0;
	double dt = T / 252;		//Time step (assuming 252 days per year)
	int num_of_steps = (int)floor(T * 252);		//Number of steps (needs to be an integer since it is used in the loop)
	double sum_of_payoffs = 0, sum_of_payoff_anti = 0;
	int flag = 0;		//Flag is used as an input in the function to return the standard normal samples, so that it generates either 4 or 2 samples per time step

	double *price1 = new double[num_of_steps + 1];		//Initializeation of price vectors 
	double *price2 = new double[num_of_steps + 1];
	double *price3 = new double[num_of_steps + 1];
	double *payoff = new double[n];						//Payoff vector of size n, 1 for each Monte Carlo simulation

	double *price1_anti = new double[num_of_steps + 1];
	double *price2_anti = new double[num_of_steps + 1];
	double *price3_anti = new double[num_of_steps + 1];
	double *payoff_anti = new double[n];

	double *payoff_average = new double[n];				//Average payoff of variables payoff and payoff_anti

	//These variables are used for the calculation of the standard deviation of the price
	double sd = 0, sd_antithetic = 0, sd_reduced = 0;

	//For the control variate method: Compare black scholes price of a european 
	//call, with the monte carlo price of the same european option and check the difference
	double price_european_BS = black_scholes(type, S0_stock1, K, T, r, vol1);
	double price_european_MC, price_paths_european, sum_of_european_payoffs = 0;
	
	double price_contol_variates; //Asian option price including the conrol variate method

	double *correl_array = new double[4];
	double fourth_sample;

	double x1, x2, x3;				//Variables for the standard normal samples
	double x1_cor, x2_cor, x3_cor;	//Correlated standard normal samples using cholesky decomposition
	double alpha;
	double x1_cor_anti, x2_cor_anti, x3_cor_anti;

	alpha = (cor23 - cor12*cor13) / sqrt(1 - cor12*cor12);

	for (j = 0; j <= n - 1; j++)
	{

		price1[0] = S0_stock1;
		price2[0] = S0_stock2;
		price3[0] = S0_stock3;

		price1_anti[0] = S0_stock1;
		price2_anti[0] = S0_stock2;
		price3_anti[0] = S0_stock3;

		for (i = 1; i <= num_of_steps; i++)
		{

			flag = (i + 1) % 2;	//Flag takes the value 0 or 1 depending on i. When flag is 0 the 
								//function twoSNormals returns 4 s normal samples other wise 4

			correl_array = twoSNormals(flag);

			if (flag == 0)
			{

				x1 = correl_array[0];
				x2 = correl_array[1];
				x3 = correl_array[2];

				fourth_sample = correl_array[3];
			}
			else
			{
				x1 = correl_array[0];
				x2 = correl_array[1];
				x3 = fourth_sample;
			}

			//Cholesky decomposition two correlate the Snormal samples
			x1_cor = x1;
			x2_cor = cor12*x1 + sqrt(1 - cor12*cor12)*x2;
			x3_cor = cor13*x1 + alpha*x2 + sqrt(1 - cor13 - alpha*alpha)*x3;

			//Antithetic samples
			x1_cor_anti = -x1_cor;
			x2_cor_anti = -x2_cor;
			x3_cor_anti = -x3_cor;

			//Price Paths using GBM model
			price1[i] = price1[i - 1] * exp((r - vol1*vol1 / 2)*dt + vol1*sqrt(dt)*x1_cor);
			price2[i] = price2[i - 1] * exp((r - vol2*vol2 / 2)*dt + vol2*sqrt(dt)*x2_cor);
			price3[i] = price3[i - 1] * exp((r - vol3*vol3 / 2)*dt + vol3*sqrt(dt)*x3_cor);

			//Antithetic simulation of the 3 stocks
			price1_anti[i] = price1_anti[i - 1] * exp((r - vol1*vol1 / 2)*dt + vol1*sqrt(dt)*x1_cor_anti);
			price2_anti[i] = price2_anti[i - 1] * exp((r - vol2*vol2 / 2)*dt + vol2*sqrt(dt)*x2_cor_anti);
			price3_anti[i] = price3_anti[i - 1] * exp((r - vol3*vol3 / 2)*dt + vol3*sqrt(dt)*x3_cor_anti);

			//Because we might not want to include in the average all the past prices I include this loop so that
			//the average variables start storing the prices for the number of days from the final day that the user
			//has requested with the past_days_average variable
			if (i > num_of_steps - past_days_average)
			{
				average1 += price1[i];
				average2 += price2[i];
				average3 += price3[i];

				//Antithetic averages
				average1_anti += price1_anti[i];
				average2_anti += price2_anti[i];
				average3_anti += price3_anti[i];
			}

		}
		
		//Calculates the average for the number of days required by the user from the final price (past_days_average)
		if (past_days_average == num_of_steps)
		{
			//If we want to average each price path over the whole period then we should include todays price which isnt sumed up during the loop
			average1 = (average1 + price1[0]) / past_days_average;
			average2 = (average2 + price2[0]) / past_days_average;
			average3 = (average3 + price3[0]) / past_days_average;

			average1_anti = (average1_anti + price1_anti[0]) / past_days_average;
			average2_anti = (average2_anti + price2_anti[0]) / past_days_average;
			average3_anti = (average3_anti + price3_anti[0]) / past_days_average;
		}
		else { //Otherwise just a simple average
			average1 = average1 / past_days_average;
			average2 = average2 / past_days_average;
			average3 = average3 / past_days_average;

			average1_anti = average1_anti / past_days_average;
			average2_anti = average2_anti / past_days_average;
			average3_anti = average3_anti / past_days_average;
		}

		//If statement for the calculation of the payoffs

		if (type == 1){		//If its a call 
			//-1 third parameter is never true, but needs to be there since the function expects 3 arguments
			payoff[j] = maxi(maxi(average1, average2, average3) - K, 0, -1);		
			payoff_anti[j] = maxi(maxi(average1_anti, average2_anti, average3_anti) - K, 0, -1);
		}
		else {			//If its a put

			payoff[j] = maxi(K - maxi(average1, average2, average3), 0, -1);		
			payoff_anti[j] = maxi(K - maxi(average1_anti, average2_anti, average3_anti), 0, -1);

		}

		//We take the average payoff from the simulation of the two payoffs for the antithetic variate variance reduction technique
		payoff_average[j] = (payoff[j] + payoff_anti[j]) / 2;

		sum_of_payoffs += payoff[j];				//Without variance reduction technique case
		sum_of_payoff_anti += payoff_average[j];	//Antithetic variate case

		average1 = 0;
		average2 = 0;
		average3 = 0;

		average1_anti = 0;
		average2_anti = 0;
		average3_anti = 0;

		
		//We also calculate a european option with similar parameters as the asian one, and thus the two will be positive
		//correlated. The difference between the european black scholes price and the european monte carlo price will be the same (or a multiple)
		//of the difference of the true value of the asian option and the value calculated by monte carlo

		// For the european option price, we only care for the final price of the underlying asset
		price_paths_european = S0_stock1 * exp((r - vol1*vol1 / 2)*T + vol1*sqrt(T)*x1);  

		if (type == 1)
			sum_of_european_payoffs += maxi(price_paths_european - K, 0, -1);
		else
			sum_of_european_payoffs += maxi(K - price_paths_european, 0, -1);

	} //End of Monte Carlo loop
	
	// Averaging the payoffs over the number of simulations and discounting back to today 
	opt_price = exp(-r*T)*sum_of_payoffs / n;
	opt_price_antithetic = exp(-r*T)*sum_of_payoff_anti / n;
	price_european_MC = exp(-r*T)*sum_of_european_payoffs / n;

	//Including the control variates variance reduction technique in the price of the asian option
	price_contol_variates = opt_price_antithetic + (price_european_BS - price_european_MC);


	//Loop for the calculation of the standard deviation of the price
	for (i = 0; i < n - 1; i++){

			sd += (payoff[i] - opt_price)*(payoff[i] - opt_price);

			sd_antithetic += (payoff_average[i] - opt_price_antithetic)*(payoff_average[i] - opt_price_antithetic);

			sd_reduced += (payoff_average[i] - price_contol_variates)*(payoff_average[i] - price_contol_variates);

	}

	//Standard Deviations for the different variance reduction techniques
	sd = sqrt(sd / (n - 1));
	sd_antithetic = sqrt(sd_antithetic / (n - 1));
	sd_reduced = sqrt(sd_reduced / (n - 1));

	cout << "Asian Option price (no variance reduction): " << opt_price << endl;
	cout << "Asian Option price with antithetic variates: " << opt_price_antithetic << endl;
	cout << endl;
	cout << "European price MC: " << price_european_MC << endl;
	cout << "European price Black scholes: " << price_european_BS << endl;
	cout << endl;
	cout << "Asian Option price with antithetic and control variates: " << price_contol_variates << endl;
	cout << endl;
	cout << "Standard Dev without variance reduction: " << sd << endl;
	cout << "Standard Dev with antithetic variates: " << sd_antithetic << endl;
	cout << "Standard Dev with antithetic and control variates: " << sd_reduced << endl;

	delete[] price1;
	delete[] price2;
	delete[] price3;
	delete[] payoff;

	delete[] price1_anti;
	delete[] price2_anti;
	delete[] price3_anti;
	delete[] payoff_anti;

	delete[] correl_array;
}

// MAIN FUNCTION
int main() {

	srand(1);

	//asianOptionPricing option(1,100, 100, 100, 100, 1, 0.05, 0.2, 0.2, 0.2, 0.75, -0.5, 0.25, 252);

	asianOptionPricing option;

	if (option.type == 1) 
		cout << "CALL OPTION!" << endl;
	else 
		cout << "PUT OPTION!" << endl;
	
	cout<< endl;
	option.priceOption();

	return 0;
}


double maxi(double a, double b, double c)
{
	double temp;

	if (a > b)
		temp = a;
	else
		temp = b;

	if (temp < c)
		temp = c;

	return temp;
}


double * twoSNormals(int flag) {


	double u1, u2, v1, v2, z1, z2, z3, z4, ss, ss2, sqrt_ss, sqrt_ss2;

	double *corr_sample = new double[4];

	double cor12 = 0.5, cor13 = -0.3, cor23 = 0.2;



	do {
		u1 = (double)rand() / RAND_MAX;
		u2 = (double)rand() / RAND_MAX;

		v1 = -1 + 2 * u1;
		v2 = -1 + 2 * u2;

		ss = v1*v1 + v2*v2;
	} while (ss > 1 || ss == 0);

	sqrt_ss = sqrt(-2 * log(ss) / ss);

	z1 = v1*sqrt_ss;
	z2 = v2*sqrt_ss;

	if (flag == 1) {
		// If flag is 1 then return only two independent standard normal samples, and use a third one from a previous draw
		corr_sample[0] = z1;
		corr_sample[1] = z2;

		return corr_sample;
	}

	do {
		u1 = (double)rand() / RAND_MAX;
		u2 = (double)rand() / RAND_MAX;

		v1 = -1 + 2 * u1;
		v2 = -1 + 2 * u2;

		ss2 = v1*v1 + v2*v2;
	} while (ss2 > 1 || ss2 == 0);

	sqrt_ss2 = sqrt(-2 * log(ss2) / ss2);

	z3 = v1*sqrt_ss2;
	z4 = v2*sqrt_ss2;


	corr_sample[0] = z1;
	corr_sample[1] = z2;
	corr_sample[2] = z3;
	corr_sample[3] = z4;

	return corr_sample;

}

// Hastings approximation to the cumulative normal
double cumulative_normal(double x){

	double const b1 = 0.31938153;
	double const b2 = -0.356563782;
	double const b3 = 1.781477937;
	double const b4 = -1.821255978;
	double const b5 = 1.330274429;
	double const p = 0.2316419;
	double const c = 0.918938533204672;

	double a, t, s, y;

	a = fabs(x);
	t = 1 / (1 + a*p);
	s = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	y = s*exp(-0.5*x*x - c);

	if (x > 0)	y = 1 - y;

	return y;

}

double black_scholes(int type_of_option,double S0, double K, double T, double r, double vol){

	double d1, d2, price;

	d1 = (log(S0 / K) + (r + vol*vol / 2)*T) / (vol*sqrt(T));

	d2 = d1 - vol*sqrt(T);

	if (type_of_option == 1 ){

		//Black Scholes Call price formula
		price = S0*cumulative_normal(d1) - K*exp(-r*T)*cumulative_normal(d2);
	}
	else {
		//Black Scholes Put price formula
		price = K*exp(-r*T)*cumulative_normal(-d2) - S0*cumulative_normal(-d1);
	}

	return price;
}
