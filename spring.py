from enum import Enum
import math

# Ported from Lua Class 'kDampedOscillator'
ZEROVEL_EPS = 0.01

class DampingType(Enum):
    STRONGLY_DAMPED = 0
    WEAKLY_DAMPED = 1
    CRITICALLY_DAMPED = 2


class DampedOscillator:
    # Create
    def __init__(self, m: float, r: float, k: float, x: float, v: float):
        self.x0_ = (x,)
        self.v0_ = v
        self.m_ = m
        self.r_ = r
        self.k_ = k
        self.xe_ = 0.0
        self.v_ = 0.0
        self.w_ = 0.0
        self.alpha = 0.0
        self.p_ = 0.0
        self.phi_ = 0.0

    def reset(self, x0: float, v0: float, xe: float):
        self.x0_ = x0 - xe
        self.v0_ = v0
        self.xe_ = xe

        # Calculate the natural angular frequency and damping ratio.

        self.w_ = math.sqrt(self.k_ / self.m_)
        self.alpha_ = self.r_ / (2 * math.sqrt(self.m_ * self.k_))

        # Calculate the solution type.

        if self.alpha_ > 1:
            self.damping_type_ = DampingType.STRONGLY_DAMPED
        elif self.alpha_ < 1:
            self.damping_type_ = DampingType.WEAKLY_DAMPED
        else:
            self.damping_type_ = DampingType.CRITICALLY_DAMPED
        # assert(false,'DT = '..self.damping_type_)
        # Solve the equation of motion.

        damping_type = self.damping_type_

        if damping_type == DampingType.STRONGLY_DAMPED:
            a = math.sqrt(self.alpha_ * self.alpha_ - 1.0)

            self.l1_ = self.w_ * (-self.alpha_ + a)
            self.l2_ = self.w_ * (-self.alpha_ - a)

            # Initial conditions.

            self.c_ = (self.v0_ - self.l1_ * self.x0_) / (self.l2_ - self.l1_)
            self.b_ = self.x0_ - self.c_

        elif damping_type == DampingType.WEAKLY_DAMPED:
            # x(t) = A*exp(-pt)*cos(vt+phi)

            self.p_ = self.r_ / self.m_ * 0.5
            self.v_ = (
                math.sqrt(4.0 * self.m_ * self.k_ - self.r_ * self.r_) / self.m_ * 0.5
            )
            self.a_ = math.sqrt(
                self.x0_ * self.x0_ + self.v0_ * self.v0_ / (self.p_ * self.p_)
            )

            if self.a_ == 0.0:
                self.phi_ = 0.0
            else:
                self.phi_ = math.acos(self.x0_ / self.a_)

        elif damping_type == DampingType.CRITICALLY_DAMPED:
            self.c_ = self.x0_
            self.b_ = self.v0_ + self.x0_ * self.w_

    def centre(self) -> float:
        return self.xe_

    def equilibrium(self, t: float, epsilon: float) -> float:
        return (
            math.abs(self.evaluate(t) - self.xe_) <= epsilon
            and math.abs(self.velocity(t)) <= epsilon
        )

    def evaluate(self, t: float) -> float:
        x = 0.0
        damping_type = self.damping_type_

        if damping_type == DampingType.STRONGLY_DAMPED:
            x = self.b_ * math.exp(self.l1_ * t) + self.c_ * math.exp(self.l2_ * t)
        elif damping_type == DampingType.WEAKLY_DAMPED:
            x = self.a_ * math.exp(-self.p_ * t) * math.cos(self.v_ * t + self.phi_)
        elif damping_type == DampingType.CRITICALLY_DAMPED:
            x = (self.b_ * t + self.c_) * math.exp(-self.w_ * t)
        else:
            print("evaluate: Exception!")

        return self.xe_ + x

    def velocity(self, t: float) -> float:
        v = 0.0
        damping_type = self.damping_type_
        if damping_type == DampingType.STRONGLY_DAMPED:
            v = self.b_ * self.l1_ * math.exp(
                self.l1_ * t
            ) + self.c_ * self.l2_ * math.exp(self.l2_ * t)
        elif damping_type == DampingType.WEAKLY_DAMPED:
            v = -self.a_ * self.p_ * math.exp(-self.p_ * t) * math.cos(
                self.v_ * t + self.phi_
            ) - self.v_ * self.a_ * math.exp(-self.p_ * t) * math.sin(
                self.v_ * t + self.phi_
            )
        elif damping_type == DampingType.CRITICALLY_DAMPED:
            v = (self.b_ - self.w_ * (self.b_ * t + self.c_)) * math.exp(-self.w_ * t)
        else:
            print("velocity: Exception!")

        return v

    def __str__(self) -> str:
        return f"{self.__class__.__name__} k:{self.k_} m:{self.m_} r:{self.r_}"



class Ktimer:

    def __init__(self, scale: float, paused: bool = False):
        self.scale_ = scale or 1.0
        self.paused_ = paused or False
        self.p1_ = 0
        self.p2_ = 0
# 
# 	reset()	-	reset the time object to now
# 
# 	The current time is used as the start time and the elapsed time is set to 0.
# 
# 	Parameters:
# 
# 		scale	-	a time scaling factor, eg 4 will slow the timer by a factor of 4
#
    def reset(self,scale: float):
        if (scale):
            self.scale_ = 1.0 / scale
        self.elapsed_time = 0.0
        self.start_time_ = self.time()

# 
# 	adjust()	-	adds a time adjustment
# 
# 	The adjustemnt is made to the start time.
# 
# 	Parameters:
# 
# 		dt	-	the time interval to add to the timer's start tiem
# 




    def adjust(self, dt: float):
        self.start_time_ = self.start_time_  + dt 

# 
# 	finished()	-	return whether the timer has finished
# 
# 	A timer is finished when the scaled elapsed time is >= 1.0.
# 
# 	Returns:
# 
# 		true if the timer has finished
# 

    def finished(self: bool):
        return self.elapsed() >= 1.0


# 
# 	range()	-	set the range for the parameterised animation value
# 
# 	The returned value is linearly interpolated and clipped to [ p1, p2 ] parametereised by
# 	the elapsed time.
# 
# 	Parameters:
# 
# 		p1, p2	-	the start and end values of the animation
#
    def  range(self,  p1: float,  p2: float ):
        self.p1_ = p1
        self.p2_ = p2 

# 
# 	p()	-	return a parameterised animation value
# 
# 	The returned value is linearly interpolated and clipped to [ p1, p2 ] parametereised by
# 	the elapsed time.
# 
# 	Parameters:
# 
# 		p1, p2	-	the start and end values of the animation
#  
# 	Returns:
# 
# 		the parameterised value
# 

    def p(self) -> float:
        return self.p2(self.p1_, self.p2_)


    def p2(self,  p1: float,  p2: float ) -> float:

        elapsed_time = self.elapsed()
        t = 0

        if (elapsed_time > 1.0):
            t = 1.0
        else:
            t = elapsed_time

        return p1 - ( p1 - p2 ) * t 
#
# 
# 	elapsed()	-	return the elapsed time
# 
# 	The elapsed time is the time that has passed since this time object was initialised or reset.
# 	Elapsed time does not include periods when the timer is paused. The elapsed time is multiplied up
# 	by the supplied scaling factor.
#
    def  elapsed(self) -> float:
    	# 
		#  If paused then return the last know elapsed time.
		#
        elapsed_time = 0

        if (self.paused_):
		 	elapsed_time = self.elapsed_time_ 
		else:
		 	elapsed_time = self.time() - self.start_time_ 
		
		return elapsed_time * self.scale_ 

# 
# 	paused()	-	accessor to the paused state of the time object
# 
# 	Returns:
# 
# 		true if the time object is paused
#
    def  paused(self) -> bool: 
		return self.paused_ 

# 
# 	pause()		-	pause the time object
#
    def  pause(self):

	# 
	#  No need to do anything if already paused.
	# 
		if ( not self.paused_ ):
			self.elapsed_time_ = self.time() - self.start_time_ 
			self.paused_ = True 

# 
# 	resume()	-	continue a paused object
# 
# 	The start time is recalculated to take into account the time
# 	at which the time object was paused.
#
    def  resume(self):
		# 
		#  No need to do anything if not paused.
		# 
		if ( self.paused_ ):
			self.start_time_ = self.time() - self.elapsed_time_ 
			self.paused_ = False


    def time(self):
        t = WorldClock.GetClock()
        t = 0.0
        return t 


class SpringBox:

    def __init__(self, startx, y, endx, m, r, k, time):
            self.osc             =   DampedOscillator( m or 1,  r or 12,  k or 100,0,0),
            self.osctimer        =   Ktimer(1,true),
            self.y               =   y,
            self.time                    =       0,
            self.finishby                =       time or 5,

    def start(self):
        self.osctimer.resume()

    def isFinished(self) -> bool:
        return (self.finishby and self.time > self.finishby)  


    def isStationary(self) -> float:
        return self.osc.equilibrium(self.osctimer.elapsed(), ZEROVEL_EPS)


    def  update(self,dt: float):
        self.time = self.time + dt

    def getPosition(self):
        elapsed = self.osctimer.elapsed()
        cx = self.osc.evaluate(elapsed)
        return cx, self.y



spring = DampedOscillator(0, 0, 0, 0, 0)

print(spring)