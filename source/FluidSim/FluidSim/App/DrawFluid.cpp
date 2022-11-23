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

	std::vector<Particle>::const_iterator itr;

	for (itr = inFluid.GetParticles().begin(); itr < inFluid.GetParticles().end(); ++itr)
	{
		mParticlePoints.push_back(itr->GetPosition());
	}
}

void DrawFluid::UpdateOpenGLBuffers()
{

}

