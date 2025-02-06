#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/detail/CalculateResiduals.hpp"
#include "Acts/EventData/detail/ParameterTraits.hpp"
#include "Acts/EventData/detail/PrintParameters.hpp"
#include "Acts/Utilities/detail/Subspace.hpp"

#include <array>
#include <cstddef>
#include <iosfwd>
#include <variant>

namespace ACTSTracking {

template<typename indices_t, std::size_t kSize>
class TMeasurement {
	static constexpr std::size_t kFullSize = Acts::detail::kParametersSize<indices_t>;

	using Subspace = Acts::detail::FixedSizeSubspace<kFullSize, kSize>;

public:
	using Scalar = Acts::ActsScalar;
	using ParametersVector = Acts::ActsVector<kSize>;
	using CovarianceMatrix = Acts::ActsSquareMatrix<kSize>;
	using FullParametersVector = Acts::ActsVector<kFullSize>;
	using ProjectionMatrix = Acts::ActsMatrix<kSize, kFullSize>;
	using ExpansionMatrix = Acts::ActsMatrix<kFullSize, kSize>;

	template<typename parameters_t, typename covariance_t>
	TMeasurement(Acts::SourceLink&& source, const std::array<indices_t, kSize> &indices,
			const Eigen::MatrixBase<parameters_t> &params,
			const Eigen::MatrixBase<covariance_t> &cov) :
			m_source(std::move(source)),
			m_subspace(indices),
			m_params(params),
			m_cov(cov) {}

	TMeasurement() = delete;
	TMeasurement(const TMeasurement&) = default;
	TMeasurement(TMeasurement&&) = default;
	~TMeasurement() = default;
	TMeasurement& operator=(const TMeasurement&) = default;
	TMeasurement& operator=(TMeasurement&&) = default;

	const Acts::SourceLink& sourceLink() const {
		return m_source;
	}

	static constexpr std::size_t size() {
		return kSize;
	}

	bool contains(indices_t i) const {
		return m_subspace.contains(i);
	}

	constexpr std::array<indices_t, kSize> indices() const {
		std::array < uint8_t, kSize > subInds = m_subspace.indices();
		std::array<indices_t, kSize> inds { };
		for (std::size_t i = 0; i < kSize; i++) {
			inds[i] = static_cast<indices_t>(subInds[i]);
		}
		return inds;
	}

	const ParametersVector& parameters() const {
		return m_params;
	}

	const CovarianceMatrix& covariance() const {
		return m_cov;
	}

	ProjectionMatrix projector() const {
		return m_subspace.template projector<Scalar>();
	}

	ExpansionMatrix expander() const {
		return m_subspace.template expander<Scalar>();
	}

	ParametersVector residuals(const FullParametersVector &reference) const {
		ParametersVector res = ParametersVector::Zero();
		Acts::detail::calculateResiduals(static_cast<indices_t>(kSize),
				m_subspace.indices(), reference, m_params, res);
		return res;
	}

	std::ostream& operator<<(std::ostream &os) const {
		Acts::detail::printMeasurement(os, static_cast<indices_t>(kSize),
				m_subspace.indices().data(), m_params.data(), m_cov.data());
		return os;
	}

private:
	Acts::SourceLink m_source;
	Subspace m_subspace;
	ParametersVector m_params;
	CovarianceMatrix m_cov;
};

template<typename parameters_t, typename covariance_t, typename indices_t,
		typename ... tail_indices_t>
auto makeMeasurement(Acts::SourceLink source,
		const Eigen::MatrixBase<parameters_t> &params,
		const Eigen::MatrixBase<covariance_t> &cov, indices_t index0,
		tail_indices_t ... tailIndices) ->
				TMeasurement<indices_t, 1u + sizeof...(tail_indices_t)>
{
	using IndexContainer = std::array<indices_t, 1u + sizeof...(tail_indices_t)>;
	return { std::move(source), IndexContainer {index0, tailIndices...}, params, cov };
}

namespace detail {

template<typename indices_t, std::size_t kN, std::size_t ... kSizes>
struct VariantMeasurementGenerator :
		VariantMeasurementGenerator<indices_t, kN - 1u, kN, kSizes...> {};

template<typename indices_t, std::size_t ... kSizes>
struct VariantMeasurementGenerator<indices_t, 0u, kSizes...> {
	using Type = std::variant<TMeasurement<indices_t, kSizes>...>;
};

}  // namespace detail

template<typename indices_t>
using VariantMeasurement = typename detail::VariantMeasurementGenerator<indices_t,
		Acts::detail::kParametersSize<indices_t>>::Type;

using BoundVariantMeasurement = VariantMeasurement<Acts::BoundIndices>;

using FreeVariantMeasurement = VariantMeasurement<Acts::FreeIndices>;

using Measurement = BoundVariantMeasurement;
using MeasurementContainer = std::vector<Measurement>;

template<typename indices_t>
std::ostream& operator<<(std::ostream &os, const VariantMeasurement<indices_t> &vm)
{
	return std::visit([&](const auto &m) {
		return (os << m);
	}, vm);
}

}  // namespace Acts
