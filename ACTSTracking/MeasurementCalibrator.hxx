#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include "SourceLink.hxx"
#include "Measurement.hxx"

namespace ACTSTracking {

class MeasurementCalibrator {
 public:
  /// Construct an invalid calibrator. Required to allow copying.
  MeasurementCalibrator() = default;
  /// Construct using a user-provided container to chose measurements from.
  MeasurementCalibrator(const MeasurementContainer& measurements)
      : m_measurements(measurements) {}

  //! Find the measurement corresponding to the source link.
  /**
   * @tparam parameters_t Track parameters type
   * @param sourceLink Input source link
   * @param parameters Input track parameters (unused)
   */
  template <typename parameters_t>
  const Measurement& operator()(const SourceLink& sourceLink,
                                const parameters_t& /* parameters */) const {
    assert((sourceLink.index() < m_measurements.size()) and
           "Source link index is outside the container bounds");
    return m_measurements[sourceLink.index()];
  }

  void calibrate(const Acts::GeometryContext& gctx,
                 const Acts::CalibrationContext& cctx,
                 const Acts::SourceLink& sourceLink,
                 Acts::VectorMultiTrajectory::TrackStateProxy trackState) const
  {
    trackState.setUncalibratedSourceLink(sourceLink);
    const auto& idxSourceLink = sourceLink.get<ACTSTracking::SourceLink>();

    assert((idxSourceLink.index() < m_measurements.size()) and
           "Source link index is outside the container bounds");

    const auto meas = std::get<1>(m_measurements[idxSourceLink.index()]);
    constexpr std::size_t kMeasurementSize = decltype(meas)::size();

    trackState.allocateCalibrated(kMeasurementSize);
    trackState.calibrated<kMeasurementSize>() = meas.parameters();
    trackState.calibratedCovariance<kMeasurementSize>() = meas.covariance();
  }

 private:
  const MeasurementContainer& m_measurements;
};

}  // namespace ACTSTracking
