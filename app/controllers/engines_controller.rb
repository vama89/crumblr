require 'net/http'
require 'date'
require 'csv'
require 'Matrix'

class EnginesController < ApplicationController

	def output
				#********************************** Input ************************************
		symbol    = %w{ba bac ge msft pep pg qcom}
		frequency = %q{m} 							#'d'=daily, 'w'=weekly, 'm'=monthly				   
		amount    = 1000
		#*****************************************************************************
		data, date, open, high, low, close, vol, adjclose, n, t = stockdata(symbol, frequency)
													# Stock data
		ret    = creturn(adjclose) 					# Stock returns
		mu     = mean(ret) 							# Mean of stock returns
		sigma  = covariance(ret) 					# Covariance of stock returns

		############# NO RISK-FREE ASSET AND SHORT SELLING ALLOWED ###################
		#************************ Minimum Variance Portfolio *************************
		mu_mv, sigma_mv, theta_mv = porfoliomv(mu, sigma)
		percenttheta_mv = (100*theta_mv).round(2)
		amounttheta_mv  = (amount*theta_mv).round(2)

		@A =  "Minimum Variance Portfolio"
		
		@B = "Portfolio Mean (%): #{(mu_mv[0,0]*100).round(2)}"
		
		@C = "Portfolio Risk (%): #{(sigma_mv[0,0]*100).round(2)}"
		 
		@D = "Portfolio Allocation for (%): "

		data = []

		for i in 0..n-1
			data << "#{symbol[i]}: #{percenttheta_mv[i,0]}"
		end

		@E = data

		
		@F = "Portfolio Allocation for ($): "

		data1 = []
		for i in 0..n-1
			data1 << "#{symbol[i]}: #{amounttheta_mv[i,0]}"
		end
		
		@G = data1

		#**************************** Tangency Portfolio *****************************
		mu_tg, sigma_tg, theta_tg = portfoliotg(mu, sigma)
		percenttheta_tg = (100*theta_tg).round(2)
		amounttheta_tg  = (amount*theta_tg).round(2)
		@H = "Tangency Portfolio"
		
		@I = "Portfolio Mean (%): #{(mu_tg[0,0]*100).round(2)}"
		
		@J = "Portfolio Risk (%): #{(sigma_tg[0,0]*100).round(2)}"
		 
		@K = "Portfolio Allocation for (%): "

		data2 = []
		for i in 0..n-1
			data2 << "#{symbol[i]}: #{percenttheta_tg[i,0]}"
		end

		@L = data2 
		
		@M = "Portfolio Allocation for ($): "

		data3 = []
		for i in 0..n-1
			data3 << "#{symbol[i]}: #{amounttheta_tg[i,0]}"
		end
		@N = data3

		#**************************** Optimal Portfolio *****************************
		g = Matrix[[2]]
		 mu_opt, sigma_opt, theta_opt = portfolioopt(g, mu, sigma)
		percenttheta_opt = (100*theta_opt).round(2)
		amounttheta_opt  = (amount*theta_opt).round(2)
		@Q = "Optimal Portfolio"
		
		@R =  "Portfolio Mean (%): #{(mu_opt[0,0]*100).round(2)}"
		
		@X = "Portfolio Risk (%): #{(sigma_opt[0,0]*100).round(2)}"
	
		@T = "Portfolio Allocation for (%): "

		data3 = []
		for i in 0..n-1
			data3 << "#{symbol[i]}: #{percenttheta_opt[i,0]}"
		end
		@U = data3

		@V =  "Portfolio Allocation for ($): "
		data4 = []
		for i in 0..n-1
			data4 << "#{symbol[i]}: #{amounttheta_opt[i,0]}"
		end
		
		@W = data4

	end

	def covariance(x)

	# Computes the covariance

		tstar = x.column_count
		covx  = x*x.transpose/(tstar-1)

	end

	def creturn(x)

	# Computes the percent change or the return

		n     = x.row_count
		t     = x.column_size
		logx  = x.map { |i| Math.log(i) }
		logx1 = logx.minor(0..n,0..t-2)
		logx2 = logx.minor(0..n,1..t-1)
		ret   = logx2 - logx1
		
		return ret

	end

	def get_historical_data(symbol, frequency)

	# Collects the raw historical stock data

		enddate   = Date.parse(Time.now.to_s)
		enddate   = enddate.strftime('%Y%m%d')
		endDate   = Date.parse(enddate)
		
		startdate = '20040101'
		startDate = Date.parse(startdate)
		
		query     = "/table.csv?s=#{symbol}&g=#{frequency}"
		query.concat("&a=#{startDate.month-1}&b=#{startDate.mday}&c=#{startDate.year}")
		query.concat("&d=#{endDate.month-1}&e=#{endDate.mday}&f=#{endDate.year.to_s}")
		
		Net::HTTP.start("itable.finance.yahoo.com", 80) { |http| 
			response = http.get(query) 
			body     = response.body.chomp
				
			return [ ] if body !~ /Date,Open,High,Low,Close,Volume,Adj Close/
				
			rows = CSV.parse(body)
				
			rows.shift
			
			return rows						
		}
		
	end 

	def Hmatrix(mu, sigma)

	# Computes elements of the H and H^(-1) matrix
	#   H = [ a b ; b c] and H^(-1) = (1/d)*[ c -b ; -b a]

		n    = mu.row_count
		bar1 = Matrix.build(n, 1) {|row, col| 1}
		a    = mu.t*sigma.inv*mu
		b    = mu.t*sigma.inv*bar1
		c    = bar1.t*sigma.inv*bar1
		d    = a*c - b**2

		return a, b, c, d
		
	end


	def mean(x)

	# Computes the mean 

		tstar = x.column_count
		sumx  = x.column(0)
		
		for i in 1..tstar-1

			sumx = x.column(i) + sumx
			
			if i == tstar-1
				meanx = sumx/tstar
			end
			
		end
		
		meanx = Matrix.column_vector(meanx)
		
		return meanx

	end

	def porfoliomv(mu, sigma)

	# Computes the minimum variance portfolio. The investor prefers a portfolio
	# with the least amount of risk and is not concern about the portfolio's
	# expected return

		a, b, c, d = Hmatrix(mu, sigma)
		n          = mu.row_count
		bar1       = Matrix.build(n, 1) {|row, col| 1}
		mu_mv      = b/c
		sigma_mv   = c.inv
		ssigma_mv  = sigma_mv**0.5
		theta_mv   = sigma.inv*bar1*sigma_mv
		
		return mu_mv, ssigma_mv, theta_mv
		
	end

	def portfolioopt(g, mu, sigma)

	# Computes the optimal portfolio. The investor prefers to maximize his
	# utility function, where the utility.
	# g is the risk aversion parameter

		a, b, c, d                = Hmatrix(mu, sigma)
		
		mu_mv, sigma_mv, theta_mv = porfoliomv(mu, sigma)
		mu_tg, sigma_tg, theta_tg = portfoliotg(mu, sigma)
			
		theta_opt                 = ((b/g)*theta_tg.t).t + ((Matrix[[1]]-b/g)*theta_mv.t).t
		mu_opt                    = mu.t*theta_opt
		sigma_opt                 = theta_opt.t*sigma*theta_opt
		ssigma_opt                = sigma_opt**(0.5)
		
		return mu_opt, sigma_opt, theta_opt

	end

	def portfoliotg(mu, sigma)

	# Computes the tangency portfolio. The investor prefers to invest in the
	# portfolio with maximum Sharpe ratio, i.e., S = Mu_EF/sigma_EF. The
	# portfolio with maximum S gives the hightes expected return per unit of
	# risk.

		a, b, c, d = Hmatrix(mu, sigma)
		mu_tg      = a/b
		sigma_tg   = (a**(0.5))/b
		theta_tg   = sigma.inv*mu*(1/b)
		
		return mu_tg, sigma_tg, theta_tg

	end


	def stockdata(symbol, frequency)

	# Organizes the historical stock data
		
		data      = [ ]						# Initialize matrices
		date      = Matrix[ ]				
		open      = Matrix[ ]
		high      = Matrix[ ]
		low       = Matrix[ ]
		close     = Matrix[ ]
		vol       = Matrix[ ]
		adjclose  = Matrix[ ]
		n         = symbol.size				# Number of stocks
		t		  = data.size				# Number of time periods
		
		(0..n-1).each do |k|
			
			data         = get_historical_data(symbol[k], frequency)	
			
			if k == 0
				t 	     = data.size		# Number of time periods
			end

			datetemp     = data[0..t].map{ |row| row[0] }
			opentemp     = data[0..t].map{ |row| row[1] }
			hightemp     = data[0..t].map{ |row| row[2] }
			lowtemp      = data[0..t].map{ |row| row[3] }
			closetemp    = data[0..t].map{ |row| row[4] }
			voltemp      = data[0..t].map{ |row| row[5] }
			adjclosetemp = data[0..t].map{ |row| row[6] }
			
			datetemp     = datetemp.map!{ |i| i.to_f }.reverse
			opentemp     = opentemp.map!{ |i| i.to_f }.reverse
			hightemp	 = hightemp.map!{ |i| i.to_f }.reverse
			closetemp    = closetemp.map!{ |i| i.to_f }.reverse
			voltemp      = voltemp.map!{ |i| i.to_f }.reverse
			adjclosetemp = adjclosetemp.map!{ |i| i.to_f }.reverse
						
			date 	     = Matrix.rows(date.to_a << datetemp)
			open 	 	 = Matrix.rows(open.to_a << opentemp)
			high 	 	 = Matrix.rows(high.to_a << hightemp)
			low 	 	 = Matrix.rows(low.to_a << lowtemp)
			close 	 	 = Matrix.rows(close.to_a << closetemp)
			vol 	 	 = Matrix.rows(vol.to_a << voltemp)
			adjclose 	 = Matrix.rows(adjclose.to_a << adjclosetemp)
			
		end
		
		return data, date, open, high, low, close, vol, adjclose, n, t
		
	end

end
