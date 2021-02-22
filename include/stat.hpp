#pragma once

template<class T>
T ComputeAvgDropOne (const std::vector<T>& v)
{
  std::size_t e = v.size();
	T sum = 0.0;
	for (std::size_t i = 1; i < e; ++i) {
	  sum += v[i];
	}
	return sum/e;
}

struct StatDropFirst
{
  StatDropFirst(std::vector<double>& v)
	{
	  Setup(v);
	}

  void Setup(std::vector<double>& v)
	{
	  // Drop the first stat
	  v.erase(v.begin());
    std::size_t n = v.size();
    auto sum = v[0];
    for (size_t i = 1; i < n; ++i)
      sum += v[i];

    mean = sum / n;

    stdev = 0;
    for (size_t i = 0; i < n; ++i)
    {
      stdev += (v[i] - mean) * (v[i] - mean);
    }
    stdev /= n;
    stdev = std::sqrt(stdev);

    minimum = *std::min_element(v.begin(), v.end());
    maximum = *std::max_element(v.begin(), v.end());
	}

  StatDropFirst(std::vector<double>& v, const std::vector<double> oth_info)
	{
    Setup(v);
		other_info = oth_info;
	}

  void Write(std::ostream& os)
  {
    os << mean    << "," << stdev   << ","
       << minimum << "," << maximum << ",";

    if (other_info.size() == 0) 
		{
		  os << std::endl;
		  return;
		}
	
		for (size_t i = 0, e = other_info.size(); i < e - 1; ++i)
		  os << other_info[i] << ",";

		os << other_info.back() << std::endl;
	}

  double mean;
  double stdev;
  double minimum;
  double maximum;
	vector<double> other_info;
};
