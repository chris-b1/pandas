
from pandas import compat
from datetime cimport *
cimport util
cimport lib
cimport numpy as np
import numpy as np
import lib
from pandas import tslib
from tslib import Timedelta, Timestamp, iNaT, NaT
from tslib import have_pytz, _get_utcoffset
from tslib cimport (
    maybe_get_tz,
    _is_utc,
    _is_tzlocal,
    _get_dst_info,
    _nat_scalar_rules,
)
cdef int64_t NPY_NAT = util.get_nat()

cdef inline int _year_add_months(pandas_datetimestruct dts, int months):
    '''new year number after shifting pandas_datetimestruct number of months'''
    return dts.year + (dts.month + months - 1) / 12

cdef inline int _month_add_months(pandas_datetimestruct dts, int months):
    '''new month number after shifting pandas_datetimestruct number of months'''
    cdef int new_month = (dts.month + months) % 12
    return 12 if new_month == 0 else new_month

cdef class DateOffset2:
    """
    cython proof of concept of DateOffset
    """

    cdef public:
        int n
        int years, months, weeks, days, hours, minutes, seconds, nanoseconds

    def __init__(self, n=0, years=0, months=0, weeks=0, days=0, hours=0,
                 minutes=0, seconds=0, nanoseconds=0):
        self.n = n
        self.years = years
        self.months = months
        self.weeks = weeks
        self.days = days
        self.hours = hours
        self.minutes = minutes
        self.seconds = seconds
        self.nanoseconds = nanoseconds

    cdef int64_t _apply(self, int64_t other):
        """ internal portion of apply, operates on int64 """
        cdef:
            pandas_datetimestruct dts
            int days_in_month

        pandas_datetime_to_datetimestruct(other, PANDAS_FR_ns, &dts)
        dts.year = _year_add_months(dts, self.months)
        dts.month = _month_add_months(dts, self.months)
        #prevent day from wrapping around month end
        days_in_month = days_per_month_table[is_leapyear(dts.year)][dts.month-1]
        dts.day = min(dts.day, days_in_month)
        # would fill in the reast of the relativedelta algo here
        # ......
        return pandas_datetimestruct_to_datetime(PANDAS_FR_ns, &dts)


    cpdef apply(self, other):
        # would have to handle other datetime types
        return Timestamp(self._apply(other.value))

    cpdef apply_index(self, dti):
        cdef:
            int64_t[:] values
            int64_t[:] out
            Py_ssize_t i, N

        values = dti.asi8
        N = len(values)
        out = np.empty(N, dtype='int64')
        for i in range(N):
            if values[i] == NPY_NAT:
                out[i] = NPY_NAT
            else:
                out[i] = self._apply(values[i])

        return dti._shallow_copy(np.asarray(out))


# ctypedef enum time_res:
#     r_min = 0
#     r_microsecond
#     r_second
#     r_minute
#     r_hour
#     r_day
#     r_month
#     r_year
#     r_max = 98
#     r_invalid = 99


# cdef conversion_factor(time_res res1, time_res res2):
#     cdef:
#         time_res min_res, max_res
#         int64_t factor

#     min_res = min(res1, res2)
#     max_res = max(res1, res2)
#     factor = 1

#     if min_res == max_res:
#         return factor

#     while min_res < max_res:
#         if min_res < r_microsecond:
#             raise "Cannot convert from less than us"
#         elif min_res == r_microsecond:
#             factor *= 1000000
#             min_res = r_second
#         elif min_res == r_second:
#             factor *= 60
#             min_res = r_minute
#         elif min_res == r_minute:
#             factor *= 60
#             min_res = r_hour
#         elif min_res == r_hour:
#             factor *= 24
#             min_res = r_day
#         else:
#             raise "Cannot convert to month or year"

#     return factor

# # Logic to generate ranges
# # -----------------------------------------------------------------------------

# cdef inline int64_t weekend_adjustment(int64_t dow, int bkwd):
#     if dow > 4:                         # sat or sun?
#         if bkwd:                        # roll back 1 or 2 days
#             return (4 - dow)
#         else:                           # roll forward 2 or 1 days
#             return (7 - dow)
#     return 0

# cdef int64_t us_in_day = conversion_factor(r_microsecond, r_day)

# cdef class _Offset:
#     """
#     Base class to generate timestamps. Set the anchor, and then move offsets
#     with next & prev. Retrieve timestamp with ts attribute.
#     """
#     cdef:
#         int64_t t, dow, biz, dayoffset
#         object start
#         _TSObject ts

#     def __cinit__(self):
#         self.t=0
#         self.dow=0
#         self.biz=0
#         self.dayoffset=0

#     cpdef anchor(self, object start=None):
#         if start is not None:
#             self.start = start
#         self.ts = convert_to_tsobject(self.start, None, None)
#         self._setup()

#     cdef _setup(self):
#         pass

#     cpdef next(self):
#         pass

#     cpdef __next__(self):
#         """wrapper around next"""
#         return self.next()

#     cpdef prev(self):
#         pass

#     cdef int64_t _ts(self):
#         """
#         Access the current timestamp value, with a possible weekday
#         adjustment.
#         """
#         cdef int64_t adj

#         if self.biz != 0:
#             adj = weekend_adjustment(self.dow, self.biz < 0)
#             return self.t + us_in_day * adj
#         else:
#             return self.t

#     cdef int64_t _get_anchor(self):
#         """
#         Retrieve an anchor relating to current offset we're on.
#         """
#         return self.t - self.dayoffset * us_in_day

#     property ts:
#         def __get__(self):
#             return self._ts()

# cdef class YearOffset(_Offset):
#     """
#     Generate annual timestamps from provided start time; apply dayoffset to
#     each timestamp. If biz > 0, we choose the next business day at each time;
#     previous if < 0.

#     Parameters
#     ----------
#     dayoffset : int
#     biz : int
#     """
#     cdef:
#         int64_t y, ly

#     def __init__(self, int64_t dayoffset=0, int64_t biz=0, object anchor=None):
#         self.dayoffset = dayoffset
#         self.biz = biz

#         if anchor is not None:
#             self.anchor(anchor)

#     cdef _setup(self):
#         cdef _TSObject ts = self.ts

#         self.t = ts.value + self.dayoffset * us_in_day
#         self.y = ts.dts.year

#         self.ly = (ts.dts.month > 2 or
#                    ts.dts.month == 2 and ts.dts.day == 29)

#         if self.biz != 0:
#             self.dow = (ts_dayofweek(ts) + self.dayoffset) % 7

#     cpdef next(self):
#         cdef int64_t days

#         days = 365 + is_leapyear(self.y + self.ly)

#         self.t += days * us_in_day
#         self.y += 1

#         if self.biz != 0:
#             self.dow = (self.dow + days) % 7

#     cpdef prev(self):
#         cdef int64_t days

#         days = 365 + is_leapyear(self.y - (1-self.ly))

#         self.t -= days * us_in_day
#         self.y -= 1

#         if self.biz != 0:
#             self.dow = (self.dow - days) % 7

# cdef class MonthOffset(_Offset):
#     """
#     Generate monthly timestamps from provided start time, and apply dayoffset
#     to each timestamp.  Stride to construct strided timestamps (eg quarterly).
#     If biz > 0, we choose the next business day at each time; previous if < 0.

#     Parameters
#     ----------
#     dayoffset : int
#     stride : int, > 0
#     biz : int
#     """
#     cdef:
#         Py_ssize_t stride, ly, m
#         int64_t y

#     def __init__(self, int64_t dayoffset=0, Py_ssize_t stride=1,
#                  int64_t biz=0, object anchor=None):
#         self.dayoffset = dayoffset
#         self.stride = stride
#         self.biz = biz

#         if stride <= 0:
#             raise ValueError("Stride must be positive")

#         if anchor is not None:
#             self.anchor(anchor)

#     cdef _setup(self):
#         cdef _TSObject ts = self.ts

#         self.t = ts.value + (self.dayoffset * us_in_day)

#         # for day counting
#         self.m  = ts.dts.month - 1
#         self.y  = ts.dts.year
#         self.ly = is_leapyear(self.y)

#         if self.biz != 0:
#             self.dow = (ts_dayofweek(ts) + self.dayoffset) % 7

#     cpdef next(self):
#         cdef:
#             int64_t tmp, days
#             Py_ssize_t j

#         days = 0
#         for j in range(0, self.stride):
#             if self.m >= 12:
#                 self.m -= 12
#                 self.y += 1
#                 self.ly = is_leapyear(self.y)
#             days += days_per_month_table[self.ly][self.m]
#             self.m += 1

#         self.t += days * us_in_day

#         if self.biz != 0:
#             self.dow = (self.dow + days) % 7

#     cpdef prev(self):
#         cdef:
#             int64_t tmp, days
#             Py_ssize_t j

#         days = 0
#         for j in range(0, self.stride):
#             self.m -= 1
#             if self.m < 0:
#                 self.m += 12
#                 self.y -= 1
#                 self.ly = is_leapyear(self.y)
#             days += days_per_month_table[self.ly][self.m]

#         self.t -= days * us_in_day

#         if self.biz != 0:
#             self.dow = (self.dow - days) % 7

# cdef class DayOfMonthOffset(_Offset):
#     """
#     Generate relative monthly timestamps from month & year of provided start
#     time. For example, fridays of the third week of each month (week=3, day=4);
#     or, thursdays of the last week of each month (week=-1, day=3).

#     Parameters
#     ----------
#     week : int
#     day : int, 0 to 6
#     """
#     cdef:
#         Py_ssize_t ly, m
#         int64_t y, day, week

#     def __init__(self, int64_t week=0, int64_t day=0, object anchor=None):
#         self.week = week
#         self.day = day

#         if self.day < 0 or self.day > 6:
#             raise ValueError("Day offset must be 0 to 6")

#         if anchor is not None:
#             self.anchor(anchor)

#     cdef _setup(self):
#         cdef _TSObject ts = self.ts

#         # rewind to beginning of month
#         self.t = ts.value - (ts.dts.day - 1) * us_in_day
#         self.dow = dayofweek(ts.dts.year, ts.dts.month, 1)

#         # for day counting
#         self.m = ts.dts.month - 1
#         self.y = ts.dts.year
#         self.ly = is_leapyear(self.y)

#     cpdef next(self):
#         cdef:
#             int64_t tmp, days

#         days = days_per_month_table[self.ly][self.m]
#         self.t += days * us_in_day
#         self.dow = (self.dow + days) % 7

#         self.m += 1
#         if self.m >= 12:
#             self.m -= 12
#             self.y += 1
#             self.ly = is_leapyear(self.y)

#     cpdef prev(self):
#         cdef:
#             int64_t tmp, days

#         days = days_per_month_table[self.ly][(self.m - 1) % 12]
#         self.t -= days * us_in_day
#         self.dow = (self.dow - days) % 7

#         self.m -= 1
#         if self.m < 0:
#             self.m += 12
#             self.y -= 1
#             self.ly = is_leapyear(self.y)

#     cdef int64_t _ts(self):
#         """
#         Overwrite default adjustment
#         """
#         cdef int64_t adj = (self.week * 7) + (self.day - self.dow) % 7
#         return self.t + us_in_day * adj

# cdef class DayOffset(_Offset):
#     """
#     Generate daily timestamps beginning with first valid time >= start time. If
#     biz != 0, we skip weekends. Stride, to construct weekly timestamps.

#     Parameters
#     ----------
#     stride : int, > 0
#     biz : boolean
#     """
#     cdef:
#         Py_ssize_t stride

#     def __init__(self, int64_t stride=1, int64_t biz=0, object anchor=None):
#         self.stride = stride
#         self.biz = biz

#         if self.stride <= 0:
#             raise ValueError("Stride must be positive")

#         if anchor is not None:
#             self.anchor(anchor)

#     cdef _setup(self):
#         cdef _TSObject ts = self.ts
#         self.t = ts.value
#         if self.biz != 0:
#             self.dow = ts_dayofweek(ts)

#     cpdef next(self):
#         self.t += (self.stride * us_in_day)
#         if self.biz != 0:
#             self.dow = (self.dow + self.stride) % 7
#             if self.dow >= 5:
#                 self.t += (7 - self.dow) * us_in_day
#                 self.dow = 0

#     cpdef prev(self):
#         self.t -= (self.stride * us_in_day)
#         if self.biz != 0:
#             self.dow = (self.dow - self.stride) % 7
#             if self.dow >= 5:
#                 self.t += (4 - self.dow) * us_in_day
#                 self.dow = 4
