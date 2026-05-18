#ifndef NVECTOR_SERIAL_DAGNOSTIC_HPP
#define NVECTOR_SERIAL_DAGNOSTIC_HPP

#include "sundials/sundials_nvector.h"

N_Vector N_VNew_SerialDiagnostic(sunindextype prognostic_length, sunindextype length, SUNContext sunctx);

N_Vector N_VMake_SerialDiagnostic(sunindextype prognostic_length, sunindextype length, sunrealtype* v_data, SUNContext sunctx);

#endif
