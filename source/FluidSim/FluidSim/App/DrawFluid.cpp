#include "DrawFluid.h"
#include "glad/glad.h"

DrawFluid::DrawFluid(const Fluid& inFluid)
{
	FromFluid(inFluid);
}

void DrawFluid::FromFluid(const Fluid& inFluid)
{
	mParticlePoints.clear();
	mParticlePoints.reserve(inFluid.GetNumParticles());

	std::vector<Particle>::const_iterator pItr;

	for (pItr = inFluid.GetParticles().begin(); pItr < inFluid.GetParticles().end(); ++pItr)
	{
		mParticlePoints.push_back(pItr->GetPosition());
	}
}

void DrawFluid::UpdateOpenGLBuffers()
{

}

