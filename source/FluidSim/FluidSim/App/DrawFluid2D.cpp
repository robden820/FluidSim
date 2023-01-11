#include "DrawFluid2D.h"
#include "glad/glad.h"

DrawFluid2D::DrawFluid2D(const Fluid2D& inFluid)
{
	FromFluid(inFluid);
}

void DrawFluid2D::FromFluid(const Fluid2D& inFluid)
{
	mParticlePoints.clear();
	mParticlePoints.reserve(inFluid.GetNumParticles());

	std::vector<Particle2D>::const_iterator pItr;

	for (pItr = inFluid.GetParticles().begin(); pItr < inFluid.GetParticles().end(); ++pItr)
	{
		mParticlePoints.push_back(pItr->GetPosition());
	}
}

void DrawFluid2D::UpdateOpenGLBuffers()
{

}

