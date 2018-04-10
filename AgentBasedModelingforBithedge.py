#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  4 20:56:01 2018

@author: bennett
"""
import numpy as np



class BaseModel(object):
    '''Base class for models to simulate markets.  Models
    will contain agents as well as asset prices as attributes
    and rules for moving the simulation forward in time.'''
    
    time = 0 # A class variable for identifying what time it is in the simulation
        
class HCModel(BaseModel):
    '''Class for simulating the HC market.'''
    
    tot_hc = 0 # A class variable for the total number of HC in the simulation
    tot_col = 0 # A class variable for the total collateral in the system

    
    def __init__(self, btc_price, btc_vol, btc_ret):
        self.agents=[] # A list of agents in the model
        self.btc_price = btc_price # The current BTC price
        self.btc_price_history = [btc_price] # Historical BTC prices
        self.btc_ret = btc_ret # The expected yearly returns of BTC
        self.btc_vol = btc_vol # The (current) volatility of BTC prices
        self.hc_price = 0 # The current HC price
        self.hc_price_history = [0]
        self.hc_bids = [] # Bid orders for HC at this time-step
        self.hc_offers = [] # Offer orders for HC at this time-step
        self.exchange_dai = 0 # Dai held by the system/exchange
        self.exchange_hc = 0 # HC held by the system/exchange
    
    def set_btc_vol(self, vol):
        '''function for setting the volatility of BTC.'''
        self.btc_vol = vol
        
    def set_btc_returns(self, exp_returns):
        '''Function for setting the expected yearly returns
        of BTC.'''
        self.btc_ret = exp_returns
    
    def btc_update_gbm(self, ret, vol, dt, price):
        '''Function for updating the BTC price according to 
        geometric Brownian motion.  Uses Euler-Maruyama method
        for comparison to more complicated update rules for which
        no closed-form solution exists.
        
        Args:
            ret - expected yearly return of BTC expressed as decimal (float)
            vol - yearly volatility of BTC expressed as decimal (float)
            dt -  time increment expressed as fraction of a year (float)
            price - current price of BTC
        '''
        
        w = np.random.normal(0, dt)
        self.btc_price = self.btc_price * (1 + ret * dt +vol *w)
        self.btc_price_history.append(self.btc_price)
        
    def debt_updates_rule1(self, btc_price_history):
        '''Function for updating the debt owed by all CDPs according to 
        our first proposed dynamic debt rule:
            D(t) = P(t0) * D(t0) / P(t) 
        where D is the debt and P is the price of BTC.
        
        Args:
            btc_price_history - [P(0), ..., P(t)] (list of floats)
        '''
        
        for agent in self.agents:
            for cdp in agent.debts:
                cdp.debt = cdp.debt * btc_price_history[-2]/btc_price_history[-1]
                
    def bid_hc(self, amt, price):
        '''Function to submit bids for HC. Requires number of HC 
        desired and at what price price.'''

    def offer_hc(self, amt, price):
        '''Function to submit offers for HC.  Requires number of
        HC offered and at what price.'''



class BaseAgent(object):
    '''Base class for agents that will participate in the market.
    Specific assets and rules of behavior will be specified in 
    the classes of particular types of agents but assets and actions
    shared by all relevant agents will be put here for convenience.'''

    def __init__(self, dai, hc, model, rule):
        self.dai = dai # Asset: Stable currency for use as collateral
        self.hc = hc # Asset: HC used to hedge BTC position
        self.model = model # The model to which this agent belongs
        self.debts = [] # A list of HC denominated debts in sublists of (HC debt, collateral, creation time)
        self.orders = [] # A list of current and previous orders
        self.rule = rule # A rule object that takes agent and model variables as inputs and outputs a (possibly empty) order
        

class HCAgent(BaseAgent):
    '''This class is for speculators who will trade HC and Dai 
    on the market in order to make a profit.  They can also print
    HC with a CDP.'''
    
    def __init__(self, dai, hc, model, rule):
        BaseAgent.__init__(self, dai, hc, model, rule)
    
    def print_hc(self, collat, hc_amt, time):
        '''Function for borrowing HC via a CDP'''
        cdp = CDP(hc_amt, collat, time)
        self.dai = self.dai - collat
        self.hc += hc_amt
        self.debts.append(cdp) 
    
    def add_collateral(self, dai):
        '''Function to add collateral in dai to a particular CDP.'''
        cdp = self.debts[debt_num]
        self.dai = self.dai - dai
        cdp.collat += dai
    
    def repay_hc(self, hc_amt):
        '''Function to repay debts, oldest first, up to a given amount
        of HC.'''
        
    def repay_debt(self, hc, debt_num):
        '''Function to repay a particular debt by its index in the debts list.'''
        cdp = self.debts[debt_num]
        self.hc = self.hc - hc
        cdp.debt = cdp.debt - hc
        if cdp.debt == 0:
            self.dai += cdp.collat # Return collateral
            del self.debts[debt_num] # Remove any cleared debts



class BTCAgent(BaseAgent):
    '''This class is for holders of BTC who will purchase HC 
    as a hedge on their BTC holdings.  Since, as a first 
    approximation, we do not implement a BTC market but 
    rather simulate the price stochastically, this agent
    will be assumed to have a static amount of BTC and will
    only have trading rules for HC. Since these agents only 
    wish to hedge their BTC positions, they don't print HC
    and won't need any CDP-related functions.'''
    
    def __init__(self, dai, btc, hc, model, rule):
        BaseAgent.__init__(self, dai, hc, model, rule)
        self.btc = btc # Asset: BTC to be hedged
 



class Order(object):
    '''Class for bid/offer orders.  Will only be for HC market at 
    first but could be extended.'''
    
    def __init__(self, otype, price, amount, time):
        self.otype = otype
        self.price = price
        self.amount = amount
        self.time = time
        self.filled = False

class CDP(object):
    '''Class for collateralized debt positions'''
    
    rate = 0.01 # A class variable for the interest rate on CDP debt (won't be used for first approximation simulations)
    lratio = 1.3 # A class variable for the liquidation ratio.  When collat < debt * lratio, we liquidate the CDP
    
    def __init__(self, debt, collateral, creation_time):
        self.debt = debt
        self.collat = collateral
        self.time = creation_time


class Rule(object):
    '''Class for Rule objects that will have data as attributes
    and decision making methods that output orders.'''
    
    def __init__(self, data, agent):
        self.data = data # A list containing all relevant data for a particular rule
        self.agent = agent # Agent object that follows this rule

    
    def example_rule(self, thr=0.1, per=2880, p=0.25, t=5760, tvol=):
        '''An example rule for illustration.  For this rule, data will be a list of
        price history of BTC and the buy/sell signal will be a slope in the recent price history.
        Details of this rule are taken from:
        
        Cocco et al., J Econ Interact Coord (2017) 12:345-365
        
        Args:
            thr - Threshold price change to trigger signal
            per - period of time over which to observe for the signal
            p - Probability of placing an order given some conditional
            t - Minimum time to wait between placing orders (time to wait for a new signal)
            
        Default values for time args chosen with units of minutes in mind:
            per defaults to 2 days
            t defaults to 4 days
            '''
        dat = self.data # Assume data is a list of BTC price history
        dat = self.agent.model.btc_price_history
        prices = dat[-per:-1] # Price over the last per time periods
        emp_vol = np.std(dat[-tvol: -1]) # Empirical volatility of BTC over last tvol time periods
        times = range(per)
        reg = np.polyfit(prices, times, 1)
        slope = per * reg[0] # Average change in price every per periods of time by linear regression
        lst_ord = self.agent.orders[-1] # The previous order placed by this agent
        ctime = self.agent.model.time # Current time in the model
        wait = ctime - lst_ord.time
        roll = np.random.binomial(1, p)
        last_price = self.agent.model.hc_price_history[-1]
        if (slope > thr) & (wait > t) & (roll == 1):
            # Some code to place a buy order (send bid to order book)
            bid = last_price * np.random.normal(1.02, emp_vol) # The following is uggested by Cocco et al.
            beta = np.random.lognormal(0.25, 0.2) 
            amt = self.agent.dai * beta / self.agent.model.hc_price
            order = Order('Bid', bid, amt, ctime)
        elif (slope < -thr) & (wait > t) & (roll == 1):
            # Some code to place a sell order (send offer to order book)
            offer = last_price / np.random.normal(1.02, emp_vol) 
            beta = np.random.lognormal(0.25, 0.2)
            amt = self.agent.hc * beta 
            order = Order('Offer', offer, amt, ctime)
        else:
            order = None
            
            
        
        




      
        
        






        