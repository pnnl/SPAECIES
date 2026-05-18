#include "nvector_serial_diagnostic.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>

namespace
{
  constexpr sunrealtype ZERO = sunrealtype(0.0);
  constexpr sunrealtype HALF = sunrealtype(0.5);
  constexpr sunrealtype ONE = sunrealtype(1.0);
  constexpr sunrealtype ONEPT5 = sunrealtype(1.5);

  struct N_VectorContent_SerialDiagnostic
  {
    sunrealtype *data;
    sunindextype prognostic_length;
    sunindextype length;
    sunbooleantype own_data;
  };

  N_VectorContent_SerialDiagnostic *content(N_Vector v)
  {
    return static_cast<N_VectorContent_SerialDiagnostic *>(v->content);
  }

  sunindextype get_prognostic_length(N_Vector v)
  {
    return content(v)->prognostic_length;
  }

  sunindextype get_length(N_Vector v)
  {
    return content(v)->length;
  }

  sunrealtype *get_data(N_Vector v)
  {
    return content(v)->data;
  }

  N_Vector create(sunrealtype *v_data, sunindextype prognostic_length,
                  sunindextype length, sunbooleantype own_data, SUNContext sunctx)
  {
    if (prognostic_length < 0 || length < 0)
    {
      return nullptr;
    }

    N_Vector v = N_VNewEmpty(sunctx);
    if (v == nullptr)
    {
      return nullptr;
    }

    v->ops->nvgetvectorid = [](N_Vector)
    {
      return SUNDIALS_NVEC_SERIAL;
    };

    v->ops->nvcloneempty = [](N_Vector w) -> N_Vector
    {
      N_Vector v = N_VNewEmpty(w->sunctx);
      if (v == nullptr)
      {
        return nullptr;
      }
      if (N_VCopyOps(w, v) != SUN_SUCCESS)
      {
        N_VFreeEmpty(v);
        return nullptr;
      }
      v->content = new N_VectorContent_SerialDiagnostic{
          nullptr, get_prognostic_length(w), get_length(w), SUNFALSE};
      return v;
    };

    v->ops->nvclone = [](N_Vector w) -> N_Vector
    {
      N_Vector v = w->ops->nvcloneempty(w);
      if (v == nullptr)
      {
        return nullptr;
      }

      const sunindextype n = get_length(w);
      if (n > 0)
      {
        content(v)->data = new sunrealtype[n];
        content(v)->own_data = SUNTRUE;
      }
      return v;
    };

    v->ops->nvdestroy = [](N_Vector v)
    {
      if (v == nullptr)
      {
        return;
      }

      if (v->content != nullptr)
      {
        if (content(v)->own_data)
        {
          delete[] get_data(v);
        }
        delete content(v);
        v->content = nullptr;
      }
      N_VFreeEmpty(v);
    };

    v->ops->nvgetarraypointer = [](N_Vector v)
    {
      return get_data(v);
    };

    v->ops->nvgetlength = [](N_Vector v)
    {
      return get_prognostic_length(v);
    };

    v->ops->nvlinearsum = [](sunrealtype a, N_Vector x, sunrealtype b,
                             N_Vector y, N_Vector z)
    {
      const sunindextype n = get_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype *yd = get_data(y);
      sunrealtype *zd = get_data(z);

      for (sunindextype i = 0; i < n; i++)
      {
        zd[i] = a * xd[i] + b * yd[i];
      }
    };

    v->ops->nvconst = [](sunrealtype c, N_Vector z)
    {
      sunrealtype *zd = get_data(z);
      const sunindextype n = get_length(z);
      for (sunindextype i = 0; i < n; i++)
      {
        zd[i] = c;
      }
    };

    v->ops->nvprod = [](N_Vector x, N_Vector y, N_Vector z)
    {
      const sunindextype n = get_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype *yd = get_data(y);
      sunrealtype *zd = get_data(z);

      for (sunindextype i = 0; i < n; i++)
      {
        zd[i] = xd[i] * yd[i];
      }
    };

    v->ops->nvdiv = [](N_Vector x, N_Vector y, N_Vector z)
    {
      const sunindextype n = get_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype *yd = get_data(y);
      sunrealtype *zd = get_data(z);

      for (sunindextype i = 0; i < n; i++)
      {
        zd[i] = xd[i] / yd[i];
      }
    };

    v->ops->nvscale = [](sunrealtype c, N_Vector x, N_Vector z)
    {
      const sunindextype n = get_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype *zd = get_data(z);

      for (sunindextype i = 0; i < n; i++)
      {
        zd[i] = c * xd[i];
      }
    };

    v->ops->nvabs = [](N_Vector x, N_Vector z)
    {
      const sunindextype n = get_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype *zd = get_data(z);

      for (sunindextype i = 0; i < n; i++)
      {
        zd[i] = std::abs(xd[i]);
      }
    };

    v->ops->nvinv = [](N_Vector x, N_Vector z)
    {
      const sunindextype n = get_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype *zd = get_data(z);

      for (sunindextype i = 0; i < n; i++)
      {
        zd[i] = ONE / xd[i];
      }
    };

    v->ops->nvaddconst = [](N_Vector x, sunrealtype b, N_Vector z)
    {
      const sunindextype n = get_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype *zd = get_data(z);

      for (sunindextype i = 0; i < n; i++)
      {
        zd[i] = xd[i] + b;
      }
    };

    v->ops->nvdotprod = [](N_Vector x, N_Vector y)
    {
      sunrealtype sum = ZERO;
      const sunindextype n = get_prognostic_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype *yd = get_data(y);

      for (sunindextype i = 0; i < n; i++)
      {
        sum += xd[i] * yd[i];
      }
      return sum;
    };

    v->ops->nvmaxnorm = [](N_Vector x)
    {
      sunrealtype max = ZERO;
      const sunindextype n = get_prognostic_length(x);
      sunrealtype *xd = get_data(x);

      for (sunindextype i = 0; i < n; i++)
      {
        max = std::max(max, std::abs(xd[i]));
      }
      return max;
    };

    v->ops->nvwrmsnorm = [](N_Vector x, N_Vector w)
    {
      sunrealtype sum = ZERO;
      const sunindextype n = get_prognostic_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype *wd = get_data(w);

      for (sunindextype i = 0; i < n; i++)
      {
        const sunrealtype prodi = xd[i] * wd[i];
        sum += prodi * prodi;
      }
      return std::sqrt(sum / n);
    };

    v->ops->nvwrmsnormmask = [](N_Vector x, N_Vector w, N_Vector id)
    {
      sunrealtype sum = ZERO;
      const sunindextype n = get_prognostic_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype *wd = get_data(w);
      sunrealtype *idd = get_data(id);

      for (sunindextype i = 0; i < n; i++)
      {
        if (idd[i] > ZERO)
        {
          const sunrealtype prodi = xd[i] * wd[i];
          sum += prodi * prodi;
        }
      }
      return std::sqrt(sum / n);
    };

    v->ops->nvmin = [](N_Vector x)
    {
      const sunindextype n = get_prognostic_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype min = xd[0];

      for (sunindextype i = 1; i < n; i++)
      {
        min = std::min(min, xd[i]);
      }
      return min;
    };

    v->ops->nvwl2norm = [](N_Vector x, N_Vector w)
    {
      sunrealtype sum = ZERO;
      const sunindextype n = get_prognostic_length(x);
      sunrealtype *xd = get_data(x);
      sunrealtype *wd = get_data(w);

      for (sunindextype i = 0; i < n; i++)
      {
        const sunrealtype prodi = xd[i] * wd[i];
        sum += prodi * prodi;
      }
      return std::sqrt(sum);
    };

    v->ops->nvl1norm = [](N_Vector x)
    {
      sunrealtype sum = ZERO;
      const sunindextype n = get_prognostic_length(x);
      sunrealtype *xd = get_data(x);

      for (sunindextype i = 0; i < n; i++)
      {
        sum += std::abs(xd[i]);
      }
      return sum;
    };

    v->ops->nvlinearcombination = [](int nvec, sunrealtype *c, N_Vector *X,
                                     N_Vector z) -> SUNErrCode
    {
      const sunindextype n = get_length(z);
      sunrealtype *zd = get_data(z);

      for (sunindextype j = 0; j < n; j++)
      {
        sunrealtype sum = ZERO;
        for (int i = 0; i < nvec; i++)
        {
          sum += c[i] * get_data(X[i])[j];
        }
        zd[j] = sum;
      }
      return SUN_SUCCESS;
    };

    v->ops->nvscaleaddmulti = [](int nvec, sunrealtype *a, N_Vector x,
                                 N_Vector *Y, N_Vector *Z) -> SUNErrCode
    {
      const sunindextype n = get_length(x);
      sunrealtype *xd = get_data(x);

      for (int i = 0; i < nvec; i++)
      {
        sunrealtype *yd = get_data(Y[i]);
        sunrealtype *zd = get_data(Z[i]);
        for (sunindextype j = 0; j < n; j++)
        {
          zd[j] = a[i] * xd[j] + yd[j];
        }
      }
      return SUN_SUCCESS;
    };

    v->ops->nvdotprodmulti = [](int nvec, N_Vector x, N_Vector *Y,
                                sunrealtype *dotprods) -> SUNErrCode
    {
      const sunindextype n = get_prognostic_length(x);
      sunrealtype *xd = get_data(x);

      for (int i = 0; i < nvec; i++)
      {
        sunrealtype *yd = get_data(Y[i]);
        dotprods[i] = ZERO;
        for (sunindextype j = 0; j < n; j++)
        {
          dotprods[i] += xd[j] * yd[j];
        }
      }
      return SUN_SUCCESS;
    };

    v->ops->nvlinearsumvectorarray = [](int nvec, sunrealtype a, N_Vector *X,
                                        sunrealtype b, N_Vector *Y,
                                        N_Vector *Z) -> SUNErrCode
    {
      const sunindextype n = get_length(Z[0]);

      for (int i = 0; i < nvec; i++)
      {
        sunrealtype *xd = get_data(X[i]);
        sunrealtype *yd = get_data(Y[i]);
        sunrealtype *zd = get_data(Z[i]);
        for (sunindextype j = 0; j < n; j++)
        {
          zd[j] = a * xd[j] + b * yd[j];
        }
      }
      return SUN_SUCCESS;
    };

    v->ops->nvscalevectorarray = [](int nvec, sunrealtype *c, N_Vector *X,
                                    N_Vector *Z) -> SUNErrCode
    {
      const sunindextype n = get_length(Z[0]);

      for (int i = 0; i < nvec; i++)
      {
        sunrealtype *xd = get_data(X[i]);
        sunrealtype *zd = get_data(Z[i]);
        for (sunindextype j = 0; j < n; j++)
        {
          zd[j] = c[i] * xd[j];
        }
      }
      return SUN_SUCCESS;
    };

    v->ops->nvconstvectorarray = [](int nvec, sunrealtype c,
                                    N_Vector *Z) -> SUNErrCode
    {
      const sunindextype n = get_length(Z[0]);

      for (int i = 0; i < nvec; i++)
      {
        sunrealtype *zd = get_data(Z[i]);
        for (sunindextype j = 0; j < n; j++)
        {
          zd[j] = c;
        }
      }
      return SUN_SUCCESS;
    };

    v->ops->nvwrmsnormvectorarray = [](int nvec, N_Vector *X, N_Vector *W,
                                       sunrealtype *nrm) -> SUNErrCode
    {
      const sunindextype n = get_length(X[0]);

      for (int i = 0; i < nvec; i++)
      {
        sunrealtype *xd = get_data(X[i]);
        sunrealtype *wd = get_data(W[i]);
        nrm[i] = ZERO;
        for (sunindextype j = 0; j < n; j++)
        {
          const sunrealtype prodi = xd[j] * wd[j];
          nrm[i] += prodi * prodi;
        }
        nrm[i] = std::sqrt(nrm[i] / n);
      }
      return SUN_SUCCESS;
    };

    v->ops->nvwrmsnormmaskvectorarray = [](int nvec, N_Vector *X, N_Vector *W,
                                           N_Vector id,
                                           sunrealtype *nrm) -> SUNErrCode
    {
      const sunindextype n = get_length(X[0]);
      sunrealtype *idd = get_data(id);

      for (int i = 0; i < nvec; i++)
      {
        sunrealtype *xd = get_data(X[i]);
        sunrealtype *wd = get_data(W[i]);
        nrm[i] = ZERO;
        for (sunindextype j = 0; j < n; j++)
        {
          if (idd[j] > ZERO)
          {
            const sunrealtype prodi = xd[j] * wd[j];
            nrm[i] += prodi * prodi;
          }
        }
        nrm[i] = std::sqrt(nrm[i] / n);
      }
      return SUN_SUCCESS;
    };

    v->ops->nvscaleaddmultivectorarray =
        [](int nvec, int nsum, sunrealtype *a, N_Vector *X, N_Vector **Y,
           N_Vector **Z) -> SUNErrCode
    {
      const sunindextype n = get_length(X[0]);

      for (int i = 0; i < nvec; i++)
      {
        sunrealtype *xd = get_data(X[i]);
        for (int j = 0; j < nsum; j++)
        {
          sunrealtype *yd = get_data(Y[j][i]);
          sunrealtype *zd = get_data(Z[j][i]);
          for (sunindextype k = 0; k < n; k++)
          {
            zd[k] = a[j] * xd[k] + yd[k];
          }
        }
      }
      return SUN_SUCCESS;
    };

    v->ops->nvlinearcombinationvectorarray =
        [](int nvec, int nsum, sunrealtype *c, N_Vector **X,
           N_Vector *Z) -> SUNErrCode
    {
      const sunindextype n = get_length(Z[0]);

      for (int j = 0; j < nvec; j++)
      {
        sunrealtype *zd = get_data(Z[j]);
        for (sunindextype k = 0; k < n; k++)
        {
          sunrealtype sum = ZERO;
          for (int i = 0; i < nsum; i++)
          {
            sum += c[i] * get_data(X[i][j])[k];
          }
          zd[k] = sum;
        }
      }
      return SUN_SUCCESS;
    };

    v->content = new N_VectorContent_SerialDiagnostic{
        v_data, prognostic_length, length, own_data};

    return v;
  }
}

N_Vector N_VNew_SerialDiagnostic(sunindextype prognostic_length, sunindextype length, SUNContext sunctx)
{
  return create(new sunrealtype[length], prognostic_length, length, SUNTRUE, sunctx);
}

N_Vector N_VMake_SerialDiagnostic(sunindextype prognostic_length, sunindextype length, sunrealtype *v_data,
                                  SUNContext sunctx)
{
  return create(v_data, prognostic_length, length, SUNFALSE, sunctx);
}
