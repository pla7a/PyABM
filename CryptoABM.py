#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 20:56:01 2018

@author: bennett
"""
import numpy as np
import math



class BaseModel(object):
    '''Base class for models to simulate markets.  Models
    will contain agents as well as asset prices as attributes
    and rules for moving the simulation forward in time.'''
    
        
class HCModel(BaseModel):
    '''Class for simulating the HC market.'''

    def __init__(self, btc_price, btc_vol, btc_ret):
        self.btc_price = btc_price # The current BTC price
        self.btc_price_history = [btc_price] # Historical BTC prices
        self.btc_ret = btc_ret # The expected yearly returns of BTC
        self.btc_vol = btc_vol # The (current) volatility of BTC prices
        self.agents = [] # A list of agents in the model
        self.new_agents = [] # Pre-generated agents to enter the simulation randomly
        self.hc_price = 1 # The current HC price
        self.hc_price_history = [1]
        self.hc_orders = [] # A list of open orders
        self.hc_bids = [] # Bid orders for HC active at this time-step. Format [(agent_id, order_id), price, amount]
        self.hc_offers = [] # Offer orders for HC active at this time-step. Format [(agent_id, order_id), price, amount]
        self.exchange_dai = 0 # Dai held by the system/exchange
        self.exchange_hc = 0 # HC held by the system/exchange
        self.exchange_cdp = [] # A list of CDPs currently under liquidation (held by exchange)
        self.liquidation_orders = [] # A list of HC bid orders related to liquidation process
        self.time = 0 # Time in the simualtion
        self.agent_count = 0 # Total number of agents ever created; used to assign agent IDs
        self.tot_hc = 0 # Total number of HC in the simulation
        self.tot_col = 0 # Total amount of collateral in the system (DAI)
        self.dt = 1 # The time step for model dynamics
 
    
    def set_btc_vol(self, vol):
        '''function for setting the volatility of BTC. This is
        a parameter used as the default value for updating the 
        price of BTC via GBM.'''
        
        self.btc_vol = vol

    
    def find_agent(self, agent_id):
        '''Return agent object associated with given agent id.'''
        ag = None
        for agent in self.agents:
            if agent.id == agent_id:
                ag = agent
                break
        return ag
    
    def set_btc_returns(self, exp_returns):
        '''Function for setting the expected yearly returns
        of BTC.'''

        self.btc_ret = exp_returns


    def agent_init_pwrlaw(self, alpha, nb, nh, maxwb, maxwh, p, rule, agent_list = 'agents'):
        '''Function to initialize the model with nb BTC agents
        and nh HC agents.  Wealth endowments are determined using
        a pair of power law distributions using exponent alpha
        and max wealth maxwb and maxwh respectively.  BTC agents
        have an average proportion p of their wealth in BTC and
        use the remainder to hedge it this exposure.
        
        Currently written with a single rule for all agents but 
        this will be changed once we've written more rules.  We
        will select from among the rules randomly in some way.'''
        
        agents = []
        for i in range(nb):
            dai = maxwb * i**(-alpha)
            frac = min(lognormal_factor(p, 0.1), 0.9)
            btc = dai * (1 - frac) / self.btc_price
            dai = dai * frac
            agent = BTCAgent(dai, btc, 0, self, rule, self.agent_count+1)
            agents.append(agent)
            self.agent_count += 1
        for i in range(nh):
            dai = maxwh * i**(-alpha)
            agent = HCAgent(dai, 0, self,rule, self.agent_count+1)
            agents.append(agent)
            self.agent_count += 1
        setattr(self, agent_list, agents)
       

    def submit_order(self, order):
        '''Function that takes an order object and submits it to the model,
        updating the order book and list of bid/offer prices accordingly.'''
        
        otype = order.otype
        if otype == 'Bid':
            self.hc_orders.append(order)
            self.hc_bids.append(order.id, order.price, order.amount)
        elif otype == 'Offer':
            self.hc_orders.append(order)
            self.hc_offers.append(order.id, order.price, order.amount)
        else:
            raise ValueError('Unrecognized order type. Neither Bid nor Offer')
    
    
    def btc_update_gbm(self, ret=None, vol=None, dt=None):
        '''Function for updating the BTC price according to 
        geometric Brownian motion.  Uses Euler-Maruyama method
        for comparison to more complicated update rules for which
        no closed-form solution exists.
        
        Args:
            ret - expected yearly return of BTC expressed as decimal (float)
            vol - yearly volatility of BTC expressed as decimal (float)
            dt -  time increment expressed as fraction of a year (float)
        '''

        # Check if the optional parameters were passed, if not then we resort to the corresponding attributes associated with the instance
        if (ret is None):
            ret = self.btc_ret
        if (vol is None):
            vol = self.btc_vol
        if (dt is None):
            dt = self.dt

        w = np.random.normal(0, dt)
        self.btc_price = self.btc_price * (1 + ret * dt +vol *w)
        self.btc_price_history.append(self.btc_price)

        
    def debt_updates_rule1(self):
        '''Function for updating the debt owed by all CDPs according to 
        our first proposed dynamic debt rule:
            D(t) = P(t0) * D(t0) / P(t) 
        where D is the debt and P is the price of BTC.
        
        Args:
            btc_price_history - [P(0), ..., P(t)] (list of floats)
        '''
        
        for agent in self.agents:
            for cdp in agent.debts:
                cdp.debt = cdp.debt * self.btc_price_history[-2]/self.btc_price_history[-1]

                
    # Every round, need to:
    #   1.) Check each CDP for liquidation condition
    #       1a.) Place HC bids to liquidate any violated CDPs
    #       1b.) Pay off CDPs and return unused collateral (if any)

    
    def liquid_check(self):
        '''Function to check all active CDPs for liquidation. For all 
        CDPs that meet liquidation condition, send HC bids to cover 
        outstanding debt (plus some cushion)'''

        debt_updates_rule1() # Update all of the debt owed by each CDP

        # Check if CDP violates liquidation ratio or if the CDP has expired (assuming 7 units of time expiration)
        for agent in self.agents:
            for cdp in agent.debts:
                if (cdp.collat < cdp.debt*CDP.lratio) or ((self.time - cdp.time) >= 7): 
                    liquidate_cdp(cdp)


    def liquidate_cdp(self, order, cdp):
        ''' Function to create match order type for the CDP liquidation process'''
        og_debt = cdp.debt
        while (cdp.debt > 0 and cdp.collat > 0):
            best_offer = max(self.hc_offers, key=lambda x: x[1]) # Fetch best current offer order
            best_price = best_offer.price # Fetch best current price
            best_amount =  best_offer.amount # Amount for the lowest offer price above
        
            # Match order type for CDP liquidation

            # If the order amount isn't large enough to finish debt/collat
            if best_amount < cdp.debt:
                if cdp.collat > best_amount*best_price:
                    cdp.collat -= best_amount*best_price
                    cdp.debt -= best_amount
                    self.hc_offers.remove(best_offer)
                    self.hc_orders.remove(best_offer)
                else:
                    cdp.debt -= cdp.collat*best_price
                    cdp.collat = 0
                    best_offer.amount -= cdp.collat*best_price

            # If order amount is larger (or equal to) than needed
            elif best_amount >= cdp.debt:
                if cdp.collat > cdp.debt*best_price:
                    cdp.collat -= cdp.debt*best_price
                    best_offer.amount -= cdp.debt
                    cdp.debt = 0

                    # If the order has been fulfilled completely then we remove it from order book
                    if (best_offer.amount == 0):
                        self.hc_offers.remove(best_offer)
                        self.hc_orders.remove(best_offer)
                else:
                    cdp.debt -= cdp.collat/best_price
                    best_offer.amount -= cdp.collat/best_price
                    cdp.collat = 0

        hc_bought = og_debt - cdp.debt
        tot_hc -= hc_bought
        cdp.agent.debts.remove(cdp) # Remove the CDP from the list of agents' debts as it has been liquidated
        
        if (cdp.collat > 0):
            cdp.agent.dai += cdp.collat # Return the remaining collateral to the agent who owned the CDP


    def price_clear(self):
        '''Function to set the HC price based on open orders in the book.
        Sort the bid/offer orders seperately.'''
        transactions = []
        self.hc_bids.sort(key=lambda x:x[2])
        self.hc_offers.sort(key=lambda x:x[2], reversed=True)
        best_bid = self.hc_bids[-1][1]
        best_offer = self.hc_offers[-1][1]
        while best_offer < best_bid:
            bid_agentid = self.hc_bids[-1][0][0]
            bid_agent = self.find_agent(bid_agentid)
            offer_agentid = self.hc_offers[-1][0][0]
            offer_agent = self.find_agent(offer_agentid)
            bid_amt = self.hc_bids[-1][2]
            offer_amt = self.hc_offers[-1[2]
            remaining = abs(bid_amt - offer_amt)
            if bid_amt < offer_amt:
                
            elif bid_amt > offer_amt:
                #stuff
            elif bid_amt = offer_amt:
                #stuff
        bid = max(self.hc_bids, key=lambda x:x[1])
        offer = max(self.hc_offers, key=lambda x:x[1])
        
    
    def update_agents(self, std):
        '''Function to update the model's list of active agents.  This
        function randomly adds or removes agents from the model.  
        Removed agents are deleted and new agents come from the 
        new_agents list. The number of agents to be added/removed is chosen
        by drawing from a normal distribution with mean at the prior 
        number of total agents.  
        
        Args: std - standard deviation of the normal dist from which
        the new agent total is drawn.'''
        current = len(self.agents)
        new = int(np.random.normal(current, std))
        change = new - current
        if change < 0:
            idxs = np.random.choice(current, abs(change))
            for idx in sorted(idxs, reversed=True):
                del(self.agents[idx])
        if change > 0:
            idxs = np.random.choice(len(self.new_agents), change)
            for idx in sorted(idxs, reversed=True):
                self.agents.append(self.new_agents[idx])
                del(self.new_agents[idx])
            

    
    def time_prop(self):
        '''Function to propagate time forward one step.'''
        self.update_agents()
        self.btc_update_gbm()
        for agent in self.agents:
            agent.rule()
        self.liquid_check()
        self.price_clear()
        self.return_collat()
    

class BaseAgent(object):
    '''Base class for agents that will participate in the market.
    Specific assets and rules of behavior will be specified in 
    the classes of particular types of agents but assets and actions
    shared by all relevant agents will be put here for convenience.'''

    def __init__(self, dai, hc, model, rule, id_num):
        self.dai = dai # Asset: Stable currency for use as collateral
        self.hc = hc # Asset: HC used to hedge BTC position
        self.model = model # The model to which this agent belongs
        self.id = id_num # Identification number for the agent
        self.rule = rule # A rule object that takes agent and model variables as inputs and outputs a (possibly empty) order
        self.orders = [] # A list of current and previous orders
        self.order_count = 0 # Number of orders placed.  Used to assign order ID number
       

class HCAgent(BaseAgent):
    '''This class is for speculators who will trade HC and Dai 
    on the market in order to make a profit.  They can also print
    HC with a CDP.  Random HC traders will be implemented here also.'''
    
    def __init__(self, dai, hc, model, rule, id_num):
        BaseAgent.__init__(self, dai, hc, model, rule, id_num)
        self.debts = [] # A list of active CDPs attached to this agent
        self.cdp_count = 0 # Number of CDPs every used by this agent.  Used to assign CDP ID number
    
    def print_hc(self, collat, hc_amt, time):
        '''Function for borrowing HC via a CDP'''
        cdp_id = (self.id_num, cdp_count)
        cdp = CDP(self, hc_amt, collat, time, cdp_id)
        self.cdp_count += 1
        self.dai = self.dai - collat
        self.hc += hc_amt
        self.debts.append(cdp) 
        self.model.hc_tot += hc_amt # Add to total number of hc in circulation
    
    def add_collateral(self, dai):
        '''Function to add collateral in DAI to a particular CDP.'''
        cdp = self.debts[debt_num]
        self.dai = self.dai - dai
        cdp.collat += dai
        
    def repay_debt(self, hc, cdp_id):
        '''Function to repay a particular debt by its index in the debts list.'''
        for cdp in self.debts:
            if cdp.id == cdp_id:
                self.hc -= hc
                cdp.debt -= hc
                self.model.hc_total -= hc # Reduce total number of hc in circulation
                if cdp.debt == 0:
                    self.dai += cdp.collat # Return collateral
                    self.debts.remove(cdp) # Remove any cleared debts
                break



class BTCAgent(BaseAgent):
    '''This class is for holders of BTC who will purchase HC 
    as a hedge on their BTC holdings.  Since, as a first 
    approximation, we do not implement a BTC market but 
    rather simulate the price stochastically, this agent
    will be assumed to have a static amount of BTC and will
    only have trading rules for HC. Since these agents only 
    wish to hedge their BTC positions, they don't print HC
    and won't need any CDP-related functions.'''
    
    def __init__(self, dai, btc, hc, model, rule, id_num):
        BaseAgent.__init__(self, dai, hc, model, rule, id_num)
        self.btc = btc # Asset: BTC to be hedged
 


class Order(object):
    '''Class for bid/offer orders.  Will only be for HC market at 
    first but could be extended.'''
    
    def __init__(self, otype, price, amount, time, id_num):
        self.otype = otype # Order types are 'Offer' or 'Bid'
        self.price = price
        self.amount = amount
        self.time = time
        self.id = id_num # Format (Agent ID, Order ID)
        self.filled = False 

class CDP(object):
    '''Class for collateralized debt positions'''
    
    rate = 0.01 # A class variable for the interest rate on CDP debt (won't be used for first approximation simulations)
    lratio = 1.3 # A class variable for the liquidation ratio.  When collat < debt * lratio, we liquidate the CDP
    
    def __init__(self, agent, debt, collateral, creation_time, id_num):
        self.agent = agent # Agent that owns this CDP
        self.debt = debt # Debt in HC to free the collateral in this CDP
        self.collat = collateral # Collateral in DAI that secures this CDP
        self.time = creation_time 
        self.id = id_num # Format (Agent_ID, CDP_ID)


class Rule(object):
    '''Class for Rule objects that will have data as attributes
    and decision making methods that output orders.'''
    
    def __init__(self, agent):
        self.agent = agent # Agent object that follows this rule


    def random_trader(self, q=0.25, p=0.1, t=5760, tvol=10):
        '''A simple rule for random trading on the HC market
        
        Args:
            q - The average proportion of wealth to trade
            p - Probability of trading
            t - Minimum time to wait between placing orders
            tvol - Time period over which to calculate empirical vol of BTC
        '''
        
        dat = self.agent.model.btc_price_history
        emp_vol = np.std(dat[-tvol:-1])
        ctime = self.agent.model.time
        lst_ord = self.agent.orders[-1]
        wait = ctime = lst_ord.time
        last_price = self.agent.model.hc_price_history[-1]
        roll = np.random.binomial(1, p)
        if (roll == 1) & (wait > t):
            roll2 = np.random.binomial(1, 0.5)
            ord_id = (self.agent.id, self.agent.order_count)
            if roll2 == 1:
                price = last_price * np.random.normal(1.02, emp_vol)
                beta = lognormal_factor(0.25, 0.2)
                amt = self.agent.dai * beta / self.agent.model.hc_price
                order = Order('Bid', price, amt, ctime, ord_id)
            elif roll2 == 0:
                price = last_price / np.random.normal(1.02, emp_vol) 
                beta = lognormal_factor(0.25, 0.2)
                amt = self.agent.hc * beta 
                order = Order('Offer', price, amt, ctime, ord_id)
                


                
        
        
        
        
    def chartist1(self, thr=0.1, per=2880, p=0.5, t=5760, tvol=10):
        '''An example rule for illustration.  For this rule, data will be a list of
        price history of BTC and the buy/sell signal will be a slope in the recent price history.
        Details of this rule are taken from:
        
        Cocco et al., J Econ Interact Coord (2017) 12:345-365
        
        Args:
            thr - Threshold price change to trigger signal
            per - Period of time over which to observe for the signal
            p - Probability of placing an order given some conditional
            t - Minimum time to wait between placing orders (time to wait for a new signal)
            tvol - Time period over which the agent calculates an empirical volatility for BTC
            
        Default values for time args chosen with units of minutes in mind:
            per defaults to 2 days
            t defaults to 4 days
            '''
            
        dat = self.agent.model.btc_price_history
        prices = dat[-per:-1] # Price over the last per time periods
        emp_vol = np.std(dat[-tvol: -1]) # Empirical volatility of BTC over last tvol time periods
        times = list(range(per))
        reg = np.polyfit(prices, times, 1)
        slope = per * reg[0] # Average change in price every per periods of time by linear regression
        lst_ord = self.agent.orders[-1] # The previous order placed by this agent
        ctime = self.agent.model.time # Current time in the model
        wait = ctime - lst_ord.time
        roll = np.random.binomial(1, p)
        last_price = self.agent.model.hc_price_history[-1]
        ord_id = (self.agent.id, self.agent.order_count)
        if (slope > thr) & (wait > t) & (roll == 1):
            # Some code to place a buy order (send bid to order book)
            bid = last_price * np.random.normal(1.02, emp_vol) # The following is suggested by Cocco et al.
            beta = lognormal_factor(0.25, 0.2)
            amt = self.agent.dai * beta / self.agent.model.hc_price
            order = Order('Bid', bid, amt, ctime, ord_id)
            self.agent.orders.append(order)
            self.agent.model.submit_order(order)
            self.agent.order_count += 1
        elif (slope < -thr) & (wait > t) & (roll == 1):
            # Some code to place a sell order (send offer to order book)
            offer = last_price / np.random.normal(1.02, emp_vol) 
            beta = lognormal_factor(0.25, 0.2)
            amt = self.agent.hc * beta 
            order = Order('Offer', offer, amt, ctime, ord_id)
            self.agent.orders.append(order)
            self.agent.model.submit_order(order)
            self.agent.order_count += 1

            
      



def lognormal_factor(mean, std):
    '''An auxilliary function to generate factors from a lognormal distribution.
    output maxes at 1.'''
    
    mu = math.log(mean / math.sqrt(1+(std / mean) ** 2))
    sig = math.sqrt(math.log(1+(std / mean) ** 2))
    beta = np.random.lognormal(mu, sig)
    beta = min(beta, 1)
    return beta
        
        
